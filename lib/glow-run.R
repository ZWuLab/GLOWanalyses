#!/usr/bin/env Rscript
# ==============================================================================
# File: lib/glow-run.R
# Purpose: The dispatcher behind the `glow` command. It runs one step of a study
#   from one YAML configuration:
#
#     glow run <step> [--config <file>] [--part <name>] [--force] [--dry-run] [--wait]
#
#   A step (data_prep, training, single_variant, snv_set) is a sequence of parts,
#   and a part runs as one unit of work or as one unit per chromosome or chunk
#   (base-config/parts.yml). The dispatcher goes through the parts in order and
#   starts the stage script of each unit. Every unit keeps a completion record, so
#   a second call of the same command skips what is finished and continues with
#   what is not. The command line names the step and the configuration and carries
#   no setting: --part, --force and --dry-run are actions.
#
#   The setting run.mode chooses where the work runs. `local` runs the units one
#   after another in this process's session. `slurm` is described in
#   lib/glow-slurm.R.
#
# Called by the launcher `glow` in the GLOWanalyses directory, which is the working
# directory. Not meant to be sourced.
#
# FUNCTIONS:
#   - glow_run_main()        Parse the command line and run the step
#   - glow_run_parse()       The command line of `glow run`
#   - glow_run_units()       The units of a part, each with its command line
#   - glow_run_script()      The stage script of a part
#   - glow_run_wanted()      Whether a part belongs to this run (run_when)
#   - glow_run_local()       Run the units of a part one after another
#

source("lib/glow-stage.R")

GLOW_STEPS <- c("data_prep", "training", "single_variant", "snv_set")

glow_run_usage <- function() c(
  "Usage: glow run <step> [--config <file>] [--part <name>] [--force] [--dry-run] [--wait]",
  "",
  "  <step>      data_prep | training | single_variant | snv_set",
  "  --config    your YAML configuration. Without it, the shipped example (example.yml).",
  "  --part      run one part of the step only (the parts are listed by --dry-run).",
  "  --force     run again what is finished, and replace results made with other settings.",
  "  --dry-run   show the units of the step and what would run, and start no unit. Two things",
  "              are written all the same: for the shipped example, its pre-built files are",
  "              put into the output folder, and with run.mode slurm the job scripts are written.",
  "  --wait      with run.mode slurm: stay until the jobs of the step have finished, submit the",
  "              units that are known only then, and stop with an error when a job failed.",
  "              In the local mode the command always stays, and the flag changes nothing.",
  "",
  "Every setting is in the configuration: base-config/base_config.yml lists them all.")

# The command line of `glow run`.
#
# args  Character vector, as commandArgs(trailingOnly = TRUE) gives it.
#
# Returns list(step, config, part, force, dry_run, wait). Stops with the usage on
# anything it does not know. The step may be spelled with a hyphen.
glow_run_parse <- function(args) {
  fail <- function(...) stop(paste(c(paste0(...), "", glow_run_usage()), collapse = "\n"), call. = FALSE)
  if (length(args) < 2L || args[1L] != "run") fail("The command is `glow run <step>`.")
  step <- gsub("-", "_", args[2L], fixed = TRUE)
  if (!step %in% GLOW_STEPS) fail("Unknown step `", args[2L], "`.")
  out <- list(step = step, config = NA_character_, part = NULL, force = FALSE, dry_run = FALSE, wait = FALSE)
  i <- 3L
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("--config", "--part") && i < length(args)) {
      out[[sub("^--", "", a)]] <- args[i + 1L]; i <- i + 2L
    } else if (a %in% c("--force", "--dry-run", "--wait")) {
      out[[gsub("-", "_", sub("^--", "", a))]] <- TRUE; i <- i + 1L
    } else fail("Unknown argument `", a, "`. A setting belongs in the configuration, not on the command line.")
  }
  out
}

# Whether a part belongs to this run.
#
# A part of the table may carry `run_when`, a list of conditions that all have to
# hold (.gs_unmet() of lib/glow-stage.R reads them, as it reads the `when` of an
# input). A part without `run_when` always belongs to the run.
#
# A part whose conditions do not hold still belongs to the run while another part of
# the run reads its output from the layout (.gs_needed_by() of lib/glow-stage.R): a
# supplied path takes a part out only when it stands in for the part for every reader.
#
# lay, parts  The layout and the table of parts. Without them only the conditions
#             of the part itself are looked at.
#
# Returns list(wanted, why): `why` says which condition does not hold.
glow_run_wanted <- function(spec, cfg, region_type, lay = NULL, parts = NULL) {
  why <- .gs_unmet(spec$run_when, cfg, region_type)
  if (nzchar(why) && !is.null(parts)) {
    run <- .gs_run_new(cfg, lay, parts, region_type); run$region_type <- region_type
    if (!is.null(.gs_needed_by(run, spec))) why <- ""
  }
  list(wanted = !nzchar(why), why = why)
}

# The units of a part, each with the arguments that select it.
#
# spec         The part, from the table.
# cfg, lay     The configuration and the layout.
# region_type  The region type of an SNV-set part, else NULL.
# parts        The table of parts.
#
# Returns a list of list(unit, args): the name of the unit (NULL for a part that is
# not split) and the selector arguments of its script. The units are those that the
# completion records expect (.gs_expected_units() of lib/glow-stage.R), so what the
# command runs and what counts as finished cannot differ. A window part takes its
# chunks from the table that the preparation wrote, so it has units only after the
# preparation has run.
glow_run_units <- function(spec, cfg, lay, region_type, parts = NULL) {
  run <- .gs_run_new(cfg, lay, parts, region_type)
  units <- .gs_expected_units(run, spec, .gs_part_folder(lay, spec, region_type))
  lapply(units, function(u) list(unit = u, args = .gs_unit_args(u)))
}

# The stage script of a part. The SNV-set scan has one script per region type
# (run-gene.R, run-window.R, run-coding.R, run-noncoding.R).
glow_run_script <- function(spec, region_type) {
  sub("{region_type}", if (is.null(region_type)) "" else region_type, spec$script, fixed = TRUE)
}

# Run the units of a part one after another.
#
# Each unit is one call of the part's stage script. A unit that is complete says so
# and ends at once. The first unit that fails stops the step: the message names it,
# and a second call of the same command continues there.
#
# Returns the number of units that were started.
glow_run_local <- function(spec, units, common_args, dry_run, status_of = NULL) {
  script <- glow_run_script(spec, status_of$region_type)
  for (u in units) {
    argv <- c(script, common_args, spec$call, u$args)
    label <- paste0(spec$step, "/", spec$part, if (!is.null(u$unit)) paste0(" ", u$unit))
    if (dry_run) {
      st <- status_of$fun(u$unit)
      cat(sprintf("  %-34s %-9s %s\n", label, st$status,
                  if (st$status %in% c("complete", "reused")) st$detail else paste0(if (nzchar(st$detail)) paste0("(", st$detail, ")  ") else "", "Rscript ", paste(argv, collapse = " "))))
      next
    }
    cat(sprintf("\n=== glow run: %s ===\n", label))
    status <- system2("Rscript", shQuote(argv))
    if (!identical(status, 0L))
      stop(sprintf("glow run: %s failed (exit status %s). Nothing after it was started. Fix the cause and run the same command again: what is finished is skipped.",
                   label, status), call. = FALSE)
  }
  length(units)
}

# Parse the command line and run the step.
glow_run_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) == 0L || args[1L] %in% c("-h", "--help", "help")) { cat(glow_run_usage(), sep = "\n"); return(invisible(0L)) }
  cmd <- glow_run_parse(args)
  config <- if (is.na(cmd$config)) "example.yml" else cmd$config
  if (!file.exists(config)) stop("Configuration not found: ", config, call. = FALSE)
  if (grepl("\\.[Rr]$", config))
    stop("`", config, "` is a configuration in the form of release 0.3.6 (a config.R). `glow run` takes one YAML ",
         "configuration: see MIGRATION.md in the GLOWanalyses directory for where each field of a config.R goes.",
         call. = FALSE)
  if (!grepl("\\.ya?ml$", config, ignore.case = TRUE))
    stop("`glow run` takes a YAML configuration (.yml): ", config, call. = FALSE)
  cfg <- GLOWpipeline::glow_absolute_paths(GLOWpipeline::read_glow_config(
    config, base = "base-config/base_config.yml", registry = "base-config/settings-registry.yml", root = "."))
  .gs_check_names(cfg)             # names that would make two units write one file, before anything starts
  lay <- GLOWpipeline::glow_layout(cfg)
  # the table of parts is read once, against the registry that the configuration was read with
  parts <- GLOWpipeline::glow_read_parts("base-config/parts.yml", attr(cfg, "glow_registry"))
  of_step <- Filter(function(p) identical(p$step, cmd$step), parts)
  if (!is.null(cmd$part)) {
    if (!cmd$part %in% vapply(of_step, `[[`, character(1), "part"))
      stop("The step `", cmd$step, "` has no part `", cmd$part, "`. Its parts, in order: ",
           paste(vapply(of_step, `[[`, character(1), "part"), collapse = ", "), ".", call. = FALSE)
    of_step <- Filter(function(p) identical(p$part, cmd$part), of_step)
  }
  # A supplied aGDS folder needs a named control folder, for a part that draws PI controls and
  # belongs to this run. A study that also supplies its PI models trains none, so its B parts
  # run without a control folder.
  for (rt in (if (cmd$step == "snv_set") cfg$snv_set$region_type else list(NULL)))
    for (spec in of_step) if (glow_run_wanted(spec, cfg, rt, lay, parts)$wanted) .gs_check_control_source(cfg, spec)
  # The example's pre-built files are copied into the layout on first use, also for a dry run,
  # so that the states shown are those of a real run.
  if (identical(cfg$dataset, "example")) .gs_seed_example(.gs_run_new(cfg, lay, parts, NULL))
  if (identical(cfg$run$records, FALSE))
    stop("`run.records` is false in this configuration, which is for the tasks of a workflow engine. ",
         "`glow run` needs the completion records: leave `run.records` at its default.", call. = FALSE)
  mode <- cfg$run$mode
  cat(sprintf("glow run %s: dataset %s, analysis %s, results under %s%s\n", cmd$step, cfg$dataset, cfg$analysis,
              lay$analysis_dir, if (cmd$dry_run) "  [dry run: nothing is started]" else ""))
  if (identical(mode, "slurm")) {
    source("lib/glow-slurm.R")
    return(invisible(glow_run_slurm(cmd, config, cfg, lay, of_step, parts)))
  }
  region_types <- if (cmd$step == "snv_set") cfg$snv_set$region_type else list(NULL)
  started <- 0L
  for (rt in region_types) {
    if (!is.null(rt)) cat(sprintf("\n--- region type: %s ---\n", rt))
    for (spec in of_step) {
      want <- glow_run_wanted(spec, cfg, rt, lay, parts)
      if (!want$wanted) { cat(sprintf("  %-34s not part of this run: %s\n", paste0(spec$step, "/", spec$part), want$why)); next }
      units <- glow_run_units(spec, cfg, lay, rt)
      if (is.null(units) || !length(units)) {
        if (cmd$dry_run) { cat(sprintf("  %-34s its units are known after the preparation has run\n", paste0(spec$step, "/", spec$part))); next }
        stop("glow run: the part `", spec$part, "` has no units to run. Its chunks come from the preparation, which has not written them.", call. = FALSE)
      }
      common <- c("--config", config, if (!is.null(rt)) c("--region-type", rt), if (cmd$force) "--force")
      status_of <- list(region_type = rt, fun = local({ spec_now <- spec; rt_now <- rt
        function(unit) glow_unit_status(cfg, lay, parts, spec_now, rt_now, unit, force = cmd$force) }))
      started <- started + glow_run_local(spec, units, common, cmd$dry_run, status_of)
    }
  }
  if (!cmd$dry_run) cat(sprintf("\nglow run %s: done. Results under %s\n", cmd$step, lay$analysis_dir))
  invisible(0L)
}

if (sys.nframe() == 0L) {
  status <- tryCatch({ glow_run_main(); 0L },
                     error = function(e) { message("Error: ", conditionMessage(e)); 1L })
  quit(save = "no", status = status)
}
