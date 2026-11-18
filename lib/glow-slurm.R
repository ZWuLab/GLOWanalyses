# ==============================================================================
# File: lib/glow-slurm.R
# Purpose: The SLURM mode of the `glow` command (run.mode: slurm). One call of
#
#     glow run <step> [--config <file>] [--part <name>] [--force] [--dry-run] [--wait]
#
#   looks at every unit of the step, submits one job for each unit that has to
#   run and is not queued, and returns. A unit that reads the output of a unit
#   that is queued or submitted in the same call is submitted with a dependency on
#   that job, so a step whose units are all known runs to its end from one call. A
#   step whose later units are known only after an earlier part (the chunks of a
#   window scan after its preparation) continues when the same command is run
#   again, or under --wait, which stays and submits the later units itself.
#
#     1. The command is safe to run again: a finished unit and a queued unit are
#        left alone, and a queued unit counts as queued whatever its record says.
#        Two calls at the same time do not submit a unit twice: a call holds a
#        reservation of the dataset while it reads the queue and submits. No job
#        submits another job.
#     2. A unit checks its inputs by their completion records when it starts, so a
#        summary never runs on an unfinished scan, whatever the scheduler reports.
#        The command decides what a unit waits for by the function that the
#        completion records use (.gs_unit_reads() of lib/glow-stage.R).
#     3. A failed unit is submitted again by a new call, together with the units
#        that depend on it. Dependent jobs are submitted with
#        --kill-on-invalid-dep=yes, so they leave the queue when a job they wait
#        for fails. Under --wait, the call stops with an error when a job it
#        follows leaves the queue without a finished unit. A unit that ran with
#        --force is finished only by the record that its own job wrote.
#     4. run.slurm.max_queued holds a call back while the user has that many tasks
#        in the queue, array tasks counted one by one, as counted when the call
#        reads the queue. It is not a hard limit across datasets: two calls for two
#        datasets that read the queue in the same moment can each fill it.
#     5. The log folder exists before submission, and a log is named by its unit
#        and job.
#     6. A job sets up R as the user's shell would: a conda environment
#        (run.slurm.conda_env) or an Rscript of the user's choice
#        (run.slurm.rscript), which the job calls by its path.
#     7. A job that starts workers asks for as many cores as it starts workers
#        (shared.n_cores), and its job script states both. A job does not run when
#        a file of the configuration changed after the submission, so the settings
#        it runs with are those its resources were asked for.
#   A scheduler that does not answer is never read as an empty queue or as a
#   submitted job: the call stops.
#
# Sourced by lib/glow-run.R, from the GLOWanalyses directory.
#
# FUNCTIONS:
#   - glow_run_slurm()         One call of `glow run` in SLURM mode
#   - glow_slurm_plan()        Every unit of the step with its state and what it waits for
#   - glow_slurm_job_name()    The job name of a unit
#   - glow_slurm_job_script()  The job script of one unit
#   - glow_slurm_resources()   The resources of a part's jobs
#   - glow_slurm_options()     run.slurm.options as arguments of sbatch
#   - glow_slurm_chain()       The checksums of the files of the configuration
#   - glow_slurm_queue()       The user's tasks in the queue
#   - glow_slurm_submit()      One call of sbatch
#

# The parts that start `shared.n_cores` forked workers inside a unit.
GLOW_WORKER_PARTS <- c("single_variant.scan", "single_variant.ld_scores", "snv_set.scan")

# The options of sbatch that the command sets itself, or that would change what one
# call of sbatch does. run.slurm.options may not carry them. --cpus-per-gpu is among
# them because sbatch refuses it beside --cpus-per-task, which every submission carries.
# --time-min is not: it lets the scheduler start a job earlier with a shorter limit, never
# a longer one than the command asks for, and that trade is the user's to make.
GLOW_SBATCH_OWNED <- list(
  long  = c("job-name", "cpus-per-task", "cpus-per-gpu", "output", "error", "chdir", "array", "dependency", "kill-on-invalid-dep",
            "parsable", "wait", "wrap", "test-only", "mem", "mem-per-cpu", "mem-per-gpu", "time", "tres-per-task"),
  short = c("J", "c", "o", "e", "D", "a", "d", "W", "t"))

# The environment variables through which sbatch takes the options above. sbatch reads
# an option from its command line first, then from the environment, then from the
# #SBATCH lines of the script. So the command gives the name, the log, the cores, the
# memory and the time of a job on the command line, and starts sbatch without the
# variables below: a job name from the environment made every job carry one name, and
# the next call could not find its jobs in the queue.
GLOW_SBATCH_ENV <- c("SBATCH_JOB_NAME", "SBATCH_OUTPUT", "SBATCH_ERROR", "SBATCH_ARRAY_INX", "SBATCH_WAIT",
                     "SBATCH_MEM_PER_NODE", "SBATCH_MEM_PER_CPU", "SBATCH_MEM_PER_GPU", "SBATCH_TIMELIMIT",
                     "SBATCH_TRES_PER_TASK", "SBATCH_CPUS_PER_GPU")

# How long squeue or sbatch may take to answer, in seconds (the R option
# glow.slurm.command_timeout, for a site with a slow scheduler and for the tests).
glow_slurm_timeout <- function() {
  t <- getOption("glow.slurm.command_timeout", 120)
  if (!is.numeric(t) || length(t) != 1L || is.na(t) || t <= 0) stop("The R option glow.slurm.command_timeout has to be a positive number of seconds.", call. = FALSE)
  t
}

# The user's tasks in the queue.
#
# Returns list(tasks, ids): the number of tasks (an array task counts as one task
# each, because squeue is asked for one line per array element), and the job ids by
# job name. Stops when squeue fails: a queue that cannot be read is not an empty
# queue.
glow_slurm_queue <- function() {
  user <- Sys.getenv("USER", unname(Sys.info()[["user"]]))
  err <- tempfile("squeue-err-"); on.exit(unlink(err, expand = FALSE), add = TRUE)
  limit <- glow_slurm_timeout()
  out <- suppressWarnings(system2("squeue", c("-h", "-r", "-u", shQuote(user), "-o", shQuote("%i %j")), stdout = TRUE, stderr = err,
                                  timeout = limit))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(as.integer(status), 0L))
    stop("squeue failed (exit status ", status, if (identical(as.integer(status), 124L)) sprintf(": no answer within %s seconds", format(limit)) else "",
         "): ", paste(c(out, if (file.exists(err)) readLines(err, warn = FALSE)), collapse = " "),
         "\nWithout the queue the command cannot tell which units are queued, so it submits nothing. ",
         "Run the same command again when the scheduler answers.", call. = FALSE)
  rows <- trimws(out); rows <- rows[nzchar(rows)]
  if (any(!grepl("^[0-9][^ ]* +[^ ]", rows)))
    stop("squeue printed a line that is not a job id and a job name: ", rows[!grepl("^[0-9][^ ]* +[^ ]", rows)][1L], call. = FALSE)
  list(tasks = length(rows), ids = split(sub(" .*$", "", rows), sub("^[^ ]+ +", "", rows)))
}

# One call of sbatch.
#
# args  The arguments, the job script last.
# what  The job name, for the message.
#
# Returns the job id. Stops when sbatch fails or prints anything but a job id: a
# submission that cannot be confirmed is not a submitted job.
glow_slurm_submit <- function(args, what) {
  err <- tempfile("sbatch-err-"); on.exit(unlink(err, expand = FALSE), add = TRUE)
  limit <- glow_slurm_timeout()
  # sbatch is started without the variables of the options that the command sets (GLOW_SBATCH_ENV)
  out <- suppressWarnings(system2("env", c(paste("-u", GLOW_SBATCH_ENV), "sbatch", shQuote(args)), stdout = TRUE, stderr = err,
                                  timeout = limit))
  status <- attr(out, "status")
  said <- paste(c(out, if (file.exists(err)) readLines(err, warn = FALSE)), collapse = " ")
  if (identical(as.integer(status), 124L))
    stop("sbatch did not answer within ", format(limit), " seconds for ", what, ". The scheduler may still have taken the job: ",
         "the next call reads the queue first and leaves the job alone if it is there. ", said, call. = FALSE)
  if (!is.null(status) && !identical(as.integer(status), 0L))
    stop("sbatch failed for ", what, " (exit status ", status, "): ", said, call. = FALSE)
  if (length(out) != 1L || !grepl("^[0-9]+(;[A-Za-z0-9_.-]+)?$", out))
    stop("sbatch did not return a job id for ", what, ": ", said, call. = FALSE)
  sub(";.*$", "", out)
}

# run.slurm.options as arguments of sbatch.
#
# The setting is one line of text, as a user would type the options after `sbatch`.
# It is split at blanks, and single or double quotes keep a value with blanks
# together (--comment="my study" is one argument). The options that the command
# sets itself are refused, in their long form (also shortened, which sbatch
# accepts), their short form, and with the value attached or apart: an option on
# the command line of sbatch overrides the line of the job script, so a second job
# name would hide the job from the next call, and a second number of cores would
# break the agreement between cores and workers. Options of the site (partition,
# account, quality of service, constraints) pass.
#
# Returns the arguments.
glow_slurm_options <- function(text) {
  if (is.null(text) || !nzchar(trimws(text))) return(character(0))
  chars <- strsplit(text, "")[[1L]]
  out <- character(0); cur <- ""; quote <- ""; started <- FALSE; i <- 1L
  while (i <= length(chars)) {
    ch <- chars[i]
    if (nzchar(quote)) {
      if (ch == quote) quote <- ""
      else if (ch == "\\" && quote == "\"" && i < length(chars) && chars[i + 1L] %in% c("\"", "\\")) { cur <- paste0(cur, chars[i + 1L]); i <- i + 1L }
      else cur <- paste0(cur, ch)
    } else if (ch %in% c("'", "\"")) { quote <- ch; started <- TRUE
    } else if (ch == "\\" && i < length(chars)) { cur <- paste0(cur, chars[i + 1L]); i <- i + 1L; started <- TRUE
    } else if (grepl("[[:space:]]", ch)) { if (started) out <- c(out, cur); cur <- ""; started <- FALSE
    } else { cur <- paste0(cur, ch); started <- TRUE }
    i <- i + 1L
  }
  if (nzchar(quote)) stop("`run.slurm.options` has a quote that is not closed: ", text, call. = FALSE)
  if (started) out <- c(out, cur)
  for (tok in out) {
    owned <- NULL
    if (startsWith(tok, "--")) {
      name <- sub("=.*$", "", substring(tok, 3L))
      hit <- GLOW_SBATCH_OWNED$long[startsWith(GLOW_SBATCH_OWNED$long, name)]
      if (nzchar(name) && length(hit)) owned <- paste0("--", hit[1L])
    } else if (startsWith(tok, "-") && nchar(tok) >= 2L) {
      # short options may be written in one word: the switches without a value first
      letters_of <- strsplit(substring(tok, 2L), "")[[1L]]
      k <- 1L
      while (k <= length(letters_of) && letters_of[k] %in% c("h", "H", "k", "O", "Q", "s", "u", "v", "V")) k <- k + 1L
      if (k <= length(letters_of) && letters_of[k] %in% GLOW_SBATCH_OWNED$short) owned <- paste0("-", letters_of[k])
    }
    if (!is.null(owned))
      stop("`run.slurm.options` carries `", tok, "` (", owned, "), which the command sets itself. The job name, the log, the ",
           "working directory, the dependencies and the cores belong to the command. Memory, time and cores of a part are set ",
           "in `run.slurm.resources`. Options of the site, such as --partition or --account, are what this setting is for.",
           call. = FALSE)
  }
  out
}

# The checksums of the files of the configuration: the files of the `extends` chain
# with the base, the registry of settings and the table of parts.
#
# cfg        The configuration as it was read.
# parts_md5  The checksum of base-config/parts.yml as it was read.
#
# Returns the lines of a file for `md5sum -c`. A job checks them before it runs. The
# checksums are those that the loader took when it read the files, not those of the
# files now: the resources of the jobs are asked for the settings that were read, and
# a file edited since then has to fail the check of the job, not pass it.
glow_slurm_chain <- function(cfg, parts_md5) {
  if (!is.character(parts_md5) || length(parts_md5) != 1L || !grepl("^[0-9a-f]{32}$", parts_md5))
    stop("glow_slurm_chain(): `parts_md5` is the checksum of the table of parts as it was read (the attribute `glow_md5` of glow_read_parts()).", call. = FALSE)
  chain <- attr(cfg, "glow_chain"); reg <- attr(cfg, "glow_registry")
  files <- c(normalizePath(chain$file), normalizePath(reg$file), normalizePath("base-config/parts.yml"))
  sums <- c(chain$md5, reg$md5, parts_md5)
  keep <- !duplicated(files)
  paste0(sums[keep], "  ", files[keep])
}

# Stop when a file of the configuration is no longer the file that was read.
#
# chain  The lines of glow_slurm_chain().
glow_slurm_chain_current <- function(chain) {
  files <- sub("^[0-9a-f]+  ", "", chain); sums <- sub("  .*$", "", chain)
  now <- unname(tools::md5sum(files))
  changed <- files[is.na(now) | now != sums]
  if (length(changed))
    stop("A file of the configuration changed while the command was running: ", paste(changed, collapse = ", "),
         ". Nothing more is submitted. Run the same command again: it reads the files as they are now.", call. = FALSE)
  invisible(TRUE)
}

# Whether the job of a unit was submitted with --force, read from its job script.
#
# A job that an earlier call submitted is followed by its name, and that call's --force
# is not known here. The job script states it (the line "# forced: yes" or "# forced:
# no" that glow_slurm_job_script() writes), and a queued unit's script is not written
# again while the unit is queued. A script that is missing or does not say is taken as
# forced, which is the stricter reading for --wait.
glow_slurm_job_forced <- function(lay, p) {
  script <- file.path(.gs_output_dir(lay, p$spec$step, p$spec$part, p$rt), "jobs", paste0(p$name, ".sh"))
  if (!file.exists(script)) return(TRUE)
  said <- grep("^# forced: (yes|no)$", readLines(script, warn = FALSE), value = TRUE)
  !identical(said, "# forced: no")
}

# The settings of the SLURM mode that only this file can judge.
glow_slurm_check_settings <- function(cfg, parts) {
  s <- cfg$run$slurm
  if (!is.null(s$max_queued) && s$max_queued < 1L)
    stop("`run.slurm.max_queued` has to be 1 or more, or null for no cap: ", s$max_queued, call. = FALSE)
  if (is.null(s$poll_seconds) || s$poll_seconds <= 0)
    stop("`run.slurm.poll_seconds` has to be a positive number of seconds.", call. = FALSE)
  known <- c("default", unique(vapply(parts, `[[`, character(1), "part")))
  unknown <- setdiff(names(s$resources), known)
  if (length(unknown))
    stop("`run.slurm.resources` names ", length(unknown), " part(s) that the workflow does not have: ", paste(unknown, collapse = ", "),
         ". The names are those of `glow run <step> --dry-run`, or `default` for every part.", call. = FALSE)
  invisible(TRUE)
}

# The resources of a part's jobs.
#
# The table of parts carries the defaults that the job templates of release 0.3.6
# used (key `slurm` of a part). `run.slurm.resources` overrides them, by part name
# or for every part under the name `default`. A part that starts workers always
# asks for shared.n_cores cores.
#
# Returns list(mem, time, cpus).
glow_slurm_resources <- function(spec, cfg) {
  res <- list(mem = "16G", time = "06:00:00", cpus = 1L)
  for (over in list(spec$slurm, cfg$run$slurm$resources[["default"]], cfg$run$slurm$resources[[spec$part]]))
    for (k in names(over)) if (!is.null(over[[k]])) res[[k]] <- over[[k]]
  if (spec$id %in% GLOW_WORKER_PARTS) res$cpus <- as.integer(cfg$shared$n_cores)
  res$cpus <- as.integer(res$cpus)
  if (!grepl("^[0-9]+[KMGT]?$", res$mem) || !grepl("^[0-9][0-9:-]*$", res$time) || is.na(res$cpus) || res$cpus < 1L)
    stop("The resources of the part `", spec$part, "` are not ones sbatch reads: mem `", res$mem, "`, time `", res$time, "`, cpus `",
         res$cpus, "`. Give mem as a number with K, M, G or T, time as minutes or [days-]hours:minutes:seconds, and cpus as a ",
         "whole number (`run.slurm.resources`).", call. = FALSE)
  res
}

# The job name of a unit.
#
# The name is what the next call looks for in the queue, so it has to belong to
# one unit only. It ends in a checksum of the real path of the unit's folder with
# the part, the region type and the unit: two studies whose names run together to
# the same text (`a-b` with `c`, and `a` with `b-c`), or that have the same names
# under two output roots, get different job names. A unit of the dataset level has
# one folder for all analyses of the dataset, so its name leaves the analysis out
# and two analyses that need it share the job.
glow_slurm_job_name <- function(cfg, lay, spec, rt, unit) {
  folder <- .gs_real_path(.gs_part_folder(lay, spec, rt))
  text <- tempfile("jobname-"); on.exit(unlink(text, expand = FALSE), add = TRUE)
  writeLines(paste(c(folder, spec$part, rt, unit), collapse = "\n"), text)
  readable <- paste(c("glow", cfg$dataset, if (spec$level != "dataset") cfg$analysis, spec$step, spec$part, rt, unit), collapse = "-")
  paste0(gsub("[^A-Za-z0-9_.-]", "_", readable), "-", substr(unname(tools::md5sum(text)), 1L, 8L))
}

# The job script of one unit.
#
# name       The job name, which is also the name of the script and of the log.
# argv       The stage script and its arguments.
# spec, cfg  The part and the configuration.
# jobs_dir, logs_dir  Where the script and the log go.
# chain      The checksums of the files of the configuration (glow_slurm_chain()),
#            or NULL for no check.
#
# Returns the path of the script. The script checks the configuration, sets up R,
# changes to the GLOWanalyses directory and runs the stage script. It starts no
# other job. Every value that comes from the configuration is passed to the shell
# as a quoted word, also inside a message, so no text of a setting is run.
glow_slurm_job_script <- function(name, argv, spec, cfg, jobs_dir, logs_dir, chain = NULL) {
  res <- glow_slurm_resources(spec, cfg)
  s <- cfg$run$slurm
  dir.create(jobs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  log <- file.path(logs_dir, paste0(name, "-%j.log"))
  check <- character(0)
  if (!is.null(chain)) {
    sums <- file.path(jobs_dir, paste0(name, ".config.md5"))
    writeLines(chain, sums)
    check <- c(
      "# The configuration has to be the one this job was submitted with: its resources were asked for those settings.",
      sprintf("md5sum --quiet -c %s || { printf 'ERROR: a file of the configuration changed after this job was submitted (see the lines above; the checksums are in %%s). The job does not run. Submit again with glow run.\\n' %s >&2; exit 1; }",
              shQuote(sums), shQuote(sums)))
  }
  setup <- if (!is.null(s$rscript)) c(
      sprintf("if [ ! -x %s ]; then printf 'ERROR: run.slurm.rscript is not an executable file: %%s\\n' %s >&2; exit 1; fi", shQuote(s$rscript), shQuote(s$rscript)),
      sprintf("RSCRIPT=%s", shQuote(s$rscript)),
      "# a stage script that starts R itself finds the same R",
      sprintf("export PATH=%s:\"$PATH\"", shQuote(dirname(s$rscript))))
    else c(
      "# conda's activation scripts read unset variables, so nounset is off around them",
      "set +u",
      "if ! command -v conda >/dev/null 2>&1; then echo \"ERROR: conda is not on PATH inside the job. Submit from a shell where conda is set up, or set run.slurm.rscript.\" >&2; exit 1; fi",
      "source \"$(conda info --base)/etc/profile.d/conda.sh\"",
      sprintf("conda activate %s || { printf 'ERROR: conda activate %%s failed. Set run.slurm.conda_env or run.slurm.rscript.\\n' %s >&2; exit 1; }", shQuote(s$conda_env), shQuote(s$conda_env)),
      "set -u",
      "command -v Rscript >/dev/null 2>&1 || { echo \"ERROR: Rscript is not on PATH after the R setup.\" >&2; exit 1; }",
      "RSCRIPT=Rscript")
  lines <- c(
    "#!/usr/bin/env bash",
    sprintf("#SBATCH --job-name=%s", name),
    # sbatch is given the log on its command line as well, which also holds for a path with a blank
    if (grepl("^[A-Za-z0-9_./%+=:,@-]+$", log)) sprintf("#SBATCH --output=%s", log),
    sprintf("#SBATCH --cpus-per-task=%d", res$cpus),
    sprintf("#SBATCH --mem=%s", res$mem),
    sprintf("#SBATCH --time=%s", res$time),
    sprintf("# Written by `glow run` (lib/glow-slurm.R) for the unit %s. It asks for %d core(s)%s.", name, res$cpus,
            if (spec$id %in% GLOW_WORKER_PARTS) sprintf(", and the unit starts shared.n_cores = %d worker(s)", as.integer(cfg$shared$n_cores)) else ""),
    sprintf("# forced: %s", if ("--force" %in% argv) "yes" else "no"),    # read back by glow_slurm_job_forced()
    "set -euo pipefail",
    check,
    s$setup,
    setup,
    sprintf("cd %s", shQuote(getwd())),
    sprintf("echo \"=== %s: started $(date) on $(hostname), Rscript $(command -v \"$RSCRIPT\") ===\"", name),
    paste("\"$RSCRIPT\"", paste(shQuote(argv), collapse = " ")),
    sprintf("echo \"=== %s: finished $(date) ===\"", name))
  path <- file.path(jobs_dir, paste0(name, ".sh"))
  writeLines(lines, path)
  path
}

# Every unit of the step with its state and what it waits for.
#
# queue        What glow_slurm_queue() returned, or NULL for no queue (a dry run
#              where squeue is not installed).
# forced_done  The units that this call has already run again under --force: a
#              later round of the same call does not force them a second time.
#
# Returns a list of units, in the order of the table. Each is a list with
#   key        part id, region type and unit
#   name       the job name
#   spec, rt, unit, argv
#   state, detail   from glow_unit_status(): complete, reused, ready, waiting, refused
#   queued     the ids of the jobs of this name in the queue
#   active     whether the unit is queued or has to run: it is not finished, it is
#              forced, or a unit it reads is active (its results will be out of date)
#   after      the keys of the units of this step whose output it reads
#   lenient    those of `after` that it reads as far as they are finished
#              (`partial` in the table): it runs after them, whatever their end
#   not_yet    the parts it reads whose units are not known yet
#   uncovered  the parts it reads that are not finished and that no unit of this
#              call produces: a part of another step, or one left out by --part
# The attribute `unknown` names the parts whose units are not known in this round.
#
# What a unit reads comes from .gs_unit_reads(), the function behind the completion
# records, so the plan and the check of a starting unit cannot differ: an input that
# the user supplies by a path is not waited for, a condition that does not hold
# drops the input, and a named unit (`unit`, `units_from`) is the only one waited
# for. A part whose units are listed in a table that an earlier part of the step
# writes (the chunks of a window scan) has no units in a round in which that part is
# active: the table on disk is the one of the run before.
glow_slurm_plan <- function(cmd, config, cfg, lay, of_step, parts, queue = NULL, forced_done = character(0)) {
  region_types <- if (cmd$step == "snv_set") cfg$snv_set$region_type else list(NULL)
  plan <- list(); unknown <- character(0); unknown_ids <- character(0); active <- character(0)
  for (rt in region_types) {
    run <- .gs_run_new(cfg, lay, parts, rt); run$region_type <- rt
    of_part <- function(id) Filter(function(p) identical(p$spec$id, id) && identical(p$rt, rt), plan)
    for (spec in of_step) {
      if (!glow_run_wanted(spec, cfg, rt, lay, parts)$wanted) next
      label <- paste0(spec$part, if (!is.null(rt)) paste0(" (", rt, ")"))
      tag <- paste(c(spec$id, rt), collapse = "|")
      kind <- if (spec$units == "region") { if (identical(rt, "window")) "chunk" else "chromosome" } else spec$units
      if (kind %in% c("chunk", "chunk_mode")) {
        writers <- unlist(lapply(spec$inputs, function(inp)
          if (identical(parts[[inp$part]]$step, cmd$step) && !nzchar(.gs_unmet(inp$when, cfg, rt))) inp$part))
        pending <- any(vapply(writers, function(id) paste(c(id, rt), collapse = "|") %in% unknown_ids ||
                                any(vapply(of_part(id), function(p) p$key %in% active, logical(1))), logical(1)))
        if (pending) { unknown <- c(unknown, label); unknown_ids <- c(unknown_ids, tag); next }
      }
      units <- glow_run_units(spec, cfg, lay, rt, parts)
      if (is.null(units) || !length(units)) { unknown <- c(unknown, label); unknown_ids <- c(unknown_ids, tag); next }
      for (u in units) {
        key <- paste(c(spec$id, rt, u$unit), collapse = "|")
        name <- glow_slurm_job_name(cfg, lay, spec, rt, u$unit)
        force <- cmd$force && !(key %in% forced_done)
        status <- glow_unit_status(cfg, lay, parts, spec, rt, u$unit, force = force, run = run)
        after <- character(0); lenient <- character(0); not_yet <- character(0); uncovered <- character(0)
        for (r in .gs_unit_reads(run, spec, u$unit)) {
          if (!is.null(r$supplied)) next                       # a file of the user, not the output of a unit
          other <- parts[[r$part]]
          here <- identical(other$step, cmd$step)
          if (here && paste(c(r$part, rt), collapse = "|") %in% unknown_ids) { not_yet <- c(not_yet, other$part); next }
          mine <- if (here) Filter(function(p) is.null(r$units) || is.null(p$unit) || p$unit %in% r$units, of_part(r$part)) else list()
          keys <- vapply(mine, `[[`, character(1), "key")
          after <- c(after, keys)
          if (r$partial) lenient <- c(lenient, keys)
          if (r$part %in% status$missing && !length(mine)) uncovered <- c(uncovered, paste0(other$step, "/", other$part))
        }
        queued <- if (is.null(queue)) NULL else queue$ids[[name]]
        is_active <- length(queued) > 0L || !status$status %in% c("complete", "reused") || any(after %in% active)
        if (is_active) active <- c(active, key)
        argv <- c(glow_run_script(spec, rt), "--config", config, if (!is.null(rt)) c("--region-type", rt), spec$call, u$args,
                  if (force) "--force")
        plan[[length(plan) + 1L]] <- list(key = key, name = name, spec = spec, rt = rt, unit = u$unit, argv = argv, forced = force,
                                          state = status$status, detail = status$detail, queued = queued, active = is_active,
                                          after = unique(after), lenient = unique(lenient), not_yet = unique(not_yet),
                                          uncovered = unique(uncovered))
      }
    }
  }
  attr(plan, "unknown") <- unknown
  plan
}

# One call of `glow run` in SLURM mode.
#
# parts  The table of parts as lib/glow-run.R read it (GLOWpipeline::glow_read_parts()).
#
# Returns 0 invisibly when every unit is finished, queued or submitted, and stops
# with a message when a unit is refused by the overwrite guard, when the scheduler
# does not answer, and under --wait when a job left the queue without a finished
# unit.
glow_run_slurm <- function(cmd, config, cfg, lay, of_step, parts) {
  if (!cmd$dry_run && (!nzchar(Sys.which("sbatch")) || !nzchar(Sys.which("squeue"))))
    stop("run.mode is slurm, and `sbatch` or `squeue` is not on PATH. Use --dry-run to write the job scripts without submitting them, or set run.mode to local.", call. = FALSE)
  # `parts` is the table as the command read it, with the checksum of the bytes it was parsed
  # from. It is not read a second time here: the units, their resources and the checksum that a
  # job checks then come from one read.
  glow_slurm_check_settings(cfg, parts)
  options <- glow_slurm_options(cfg$run$slurm$options)
  chain <- glow_slurm_chain(cfg, attr(parts, "glow_md5"))
  cap <- cfg$run$slurm$max_queued
  forced_done <- character(0)
  label_of <- function(p) paste(c(p$spec$part, p$rt, p$unit), collapse = " ")
  repeat {
    # One call at a time reads the queue and submits for a dataset, so that two calls do not
    # each find a unit missing from the queue and submit it. The reservation is held for the
    # round only, not while --wait sleeps. A call waits up to two minutes for the reservation
    # of another call (the R option glow.slurm.reserve_wait, in seconds, is for the tests).
    reservation <- if (!cmd$dry_run) GLOWpipeline::glow_lock(lay$dataset_dir, "slurm-submission",
                                                             wait = getOption("glow.slurm.reserve_wait", 120))
    round <- tryCatch({
      glow_slurm_chain_current(chain)                  # the files are still those that were read
      queue <- if (cmd$dry_run && !nzchar(Sys.which("squeue"))) NULL else glow_slurm_queue()
      plan <- glow_slurm_plan(cmd, config, cfg, lay, of_step, parts, queue, forced_done)
      refused <- Filter(function(p) p$state == "refused" && !length(p$queued), plan)
      if (length(refused))
        stop("glow run: ", length(refused), " unit(s) hold finished results made with other dependencies:\n",
             paste0("  - ", vapply(refused, `[[`, character(1), "name"), ": ", vapply(refused, `[[`, character(1), "detail"), collapse = "\n"),
             "\nGive the run a new `analysis` name, or run again with --force. Nothing was submitted.", call. = FALSE)
      active <- vapply(Filter(function(p) p$active, plan), `[[`, character(1), "key")
      room <- if (is.null(cap)) Inf else max(0, cap - (if (is.null(queue)) 0L else queue$tasks))
      job_of <- list(); held <- 0L; waiting <- 0L; submitted <- 0L
      for (p in plan) {
        runs_first <- intersect(p$after, active)         # the units it reads that are queued or have to run
        if (length(p$queued)) {
          job_of[[p$key]] <- p$queued
          action <- paste("in the queue as job", paste(p$queued, collapse = ", "))
        } else if (!p$active) {
          action <- p$state
        } else if (length(p$not_yet)) {
          waiting <- waiting + 1L; action <- paste("waits for the units of", paste(p$not_yet, collapse = ", "))
        } else if (length(p$uncovered)) {
          waiting <- waiting + 1L
          action <- paste0("waits for ", paste(p$uncovered, collapse = ", "), ", which this call does not run. Run that part first.")
        } else if (!all(runs_first %in% names(job_of))) {
          waiting <- waiting + 1L; action <- "waits for a unit of this step that is not submitted yet"
        } else if (p$state == "waiting" && !length(runs_first)) {
          waiting <- waiting + 1L; action <- paste("waits", p$detail)
        } else if (room <= 0) {
          held <- held + 1L; action <- "held back by run.slurm.max_queued"
        } else {
          step_dir <- .gs_output_dir(lay, p$spec$step, p$spec$part, p$rt)
          logs_dir <- file.path(step_dir, "logs", "slurm")
          script <- glow_slurm_job_script(p$name, p$argv, p$spec, cfg, file.path(step_dir, "jobs"), logs_dir, chain)
          ok <- unlist(job_of[setdiff(runs_first, p$lenient)]); any_end <- unlist(job_of[intersect(runs_first, p$lenient)])
          dependency <- paste(c(if (length(ok)) paste0("afterok:", paste(ok, collapse = ":")),
                                if (length(any_end)) paste0("afterany:", paste(any_end, collapse = ":"))), collapse = ",")
          # what the command sets goes on the command line of sbatch, which wins over the environment
          res <- glow_slurm_resources(p$spec, cfg)
          args <- c("--parsable", paste0("--job-name=", p$name), paste0("--output=", file.path(logs_dir, paste0(p$name, "-%j.log"))),
                    sprintf("--cpus-per-task=%d", res$cpus), paste0("--mem=", res$mem), paste0("--time=", res$time),
                    if (nzchar(dependency)) c("--kill-on-invalid-dep=yes", paste0("--dependency=", dependency)), options, script)
          room <- room - 1
          if (cmd$dry_run) {
            job_of[[p$key]] <- paste0("<", p$name, ">")
            action <- paste("would submit: sbatch", paste(shQuote(args), collapse = " "))
          } else {
            id <- tryCatch(glow_slurm_submit(args, p$name), error = function(e)
              stop(conditionMessage(e), "\n", submitted, " job(s) were submitted before it and stay in the queue. ",
                   "Run the same command again when the scheduler answers: what is queued is left alone.", call. = FALSE))
            job_of[[p$key]] <- id; submitted <- submitted + 1L
            action <- paste0("submitted as job ", id, if (nzchar(dependency)) paste0(", ", dependency))
          }
        }
        cat(sprintf("  %-44s %s\n", label_of(p), action))
      }
      list(plan = plan, job_of = job_of, held = held, waiting = waiting, submitted = submitted, active = active)
    }, finally = GLOWpipeline::glow_unit_release(reservation))
    plan <- round$plan; job_of <- round$job_of
    unknown <- attr(plan, "unknown")
    if (length(unknown)) cat(sprintf("  %-44s its units are known after the part before it has run\n", unknown), sep = "")
    if (!length(round$active) && !length(unknown)) { cat("\nglow run ", cmd$step, ": every unit is finished.\n", sep = ""); return(invisible(0L)) }
    if (cmd$dry_run) { cat("\n[dry run: the job scripts are written, nothing is submitted]\n"); return(invisible(0L)) }
    if (!cmd$wait) {
      more <- length(unknown) > 0L || round$held > 0L || round$waiting > 0L
      cat(sprintf("\nglow run %s: %d job(s) submitted, %d already in the queue.%s%s\n", cmd$step, round$submitted,
                  length(job_of) - round$submitted,
                  if (more) " Run the same command again when they have finished, or use --wait, to submit what is left." else
                    " Run the same command again to see the state.",
                  if (cmd$force && length(unknown)) paste0(" The units of ", paste(unknown, collapse = ", "),
                    " are known only then. To run them again as well, use --force together with --wait, or give each of those parts",
                    " to --part with --force afterwards.") else ""))
      return(invisible(0L))
    }
    # --wait: stay until the jobs of this round have left the queue, then judge them and look again
    mine <- unlist(job_of)
    if (!length(mine)) {
      if (round$held == 0L)
        stop("glow run --wait: nothing is queued and ", round$waiting, " unit(s) cannot be submitted. See the lines above.", call. = FALSE)
      cat(sprintf("\n[--wait] the queue is full (run.slurm.max_queued = %s); looking again in %s second(s) ...\n", cap, cfg$run$slurm$poll_seconds))
      Sys.sleep(cfg$run$slurm$poll_seconds)
      next
    }
    cat(sprintf("\n[--wait] %d job(s) in the queue; looking again when they have left it ...\n", length(mine)))
    told <- FALSE
    repeat {
      Sys.sleep(cfg$run$slurm$poll_seconds)
      # a scheduler that does not answer says nothing about the jobs: keep waiting
      now <- tryCatch(glow_slurm_queue(), error = function(e) {
        if (!told) message("[--wait] ", conditionMessage(e), "\n[--wait] still waiting; the queue is read again every ", cfg$run$slurm$poll_seconds, " second(s).")
        told <<- TRUE; NULL })
      if (!is.null(now) && !any(mine %in% unlist(now$ids))) break
    }
    followed <- Filter(function(p) !is.null(job_of[[p$key]]), plan)
    # A followed unit is finished when it is complete. A unit whose job ran with --force has to
    # be complete by the record that this job wrote: the record names its job (`producer`). The
    # valid record of an earlier run is not that (a forced job that failed before it began
    # leaves it in place), nor is the record of another run, nor a folder that an analysis up
    # the chain finished, which a forced unit does not reuse. A record that names the job counts
    # whenever it was written, also before this call looked at the queue. A unit that was not
    # forced is finished by any valid record: its job skips a unit that another run completed.
    unfinished <- function(p) {
      status <- glow_unit_status(cfg, lay, parts, p$spec, p$rt, p$unit, force = FALSE)$status
      forced <- if (length(p$queued)) glow_slurm_job_forced(lay, p) else isTRUE(p$forced)
      if (!forced) return(!status %in% c("complete", "reused"))
      if (!identical(status, "complete")) return(TRUE)
      record <- GLOWpipeline::glow_unit_record(.gs_part_folder(lay, p$spec, p$rt), p$spec$part, p$unit)
      !isTRUE(unname(record$producer["slurm_job"]) %in% as.character(job_of[[p$key]]))
    }
    failed <- Filter(unfinished, followed)
    if (length(failed))
      stop("glow run --wait: ", length(failed), " of ", length(followed), " unit(s) left the queue without finishing: ",
           paste(utils::head(vapply(failed, label_of, character(1)), 10L), collapse = ", "), if (length(failed) > 10L) ", ..." else "",
           ". Their logs are under logs/slurm/ of the step folder. A unit that waited for a failed one was taken out of the queue. ",
           "Fix the cause and run the same command again: what is finished is skipped.", call. = FALSE)
    forced_done <- c(forced_done, names(job_of))       # a forced unit has run once; what follows continues from the records
  }
}
