# ==============================================================================
# File: lib/glow-stage.R
# Purpose: The common start of every GLOWanalyses stage script. One function,
#   glow_stage_begin(), parses the command line, reads the configuration and tells the
#   script where its outputs go.
#
#   The configuration is a YAML file (or none, which means example.yml). The settings are
#   read and checked by GLOWpipeline::read_glow_config(), the folders follow from the
#   three names output_root / dataset / analysis (GLOWpipeline::glow_layout()), the
#   inputs that an earlier part produces are filled from that layout, and the settings
#   are assigned into the script's environment under the variable names the scripts
#   have always read. A config.R, the form of release 0.3.6, is refused with a pointer
#   to the migration guide.
#
# Sourced (not run) by the scripts of 00-data-prep/ to 04-summary/, from the
# GLOWanalyses directory.
#
#   The helper also keeps the records of the run
#   (base-config/parts.yml says what each part depends on): a unit that is complete
#   is skipped, a run that would overwrite results made with other settings is
#   refused, and a part that the analysis extends from another analysis is reused.
#
# MAIN FUNCTIONS:
#   - glow_stage_begin()            The start of a stage script (see above)
#   - glow_stage_end()              The end of a stage script: the completion record
#   - glow_stage_option()           An option of release 0.3.6 that is a setting
#   - glow_stage_input_dir()        The finished results of this run that a gathering part reads
#   - glow_stage_watch(), glow_stage_written()   The files a part wrote into a shared folder
#   - glow_stage_clear_dir()        An empty folder for a part that owns one
#   - glow_scan_unit_files()        The files of one unit of an SNV-set scan
#   - glow_unit_status()            The state of a unit, without running it
#   - glow_link_dataset_assets()    Links from an analysis to the dataset's data
#
# INTERNAL HELPERS:
#   - .gs_parse_args()              The command line of a stage script
#   - .gs_unit_reads(), .gs_input_units()   What a unit reads, for the records and the SLURM mode
#   - .gs_needed_by()               The part of a run that still reads a part that is left out
#   - .gs_real_path()               The real path of a place that may not exist yet
#   - .gs_refuse_option()           Stop for a setting on the command line of a YAML run
#   - .gs_check_names()             Names that would make outputs collide
#   - .gs_region_type()             The one region type of an SNV-set run
#   - .gs_derived_rules()           Where each derived input comes from
#   - .gs_derive_inputs()           Inputs left null, filled from the layout
#   - .gs_effective_value()         A setting's value as the run uses it
#   - .gs_control_format(), .gs_control_source()   The control files of the PI training
#   - .gs_as_predicates()           Annotation tests as GLOWr objects
#   - .gs_assign_settings()         The settings under the scripts' variable names
#   - .gs_dataprep_paths()          The folders of the data-preparation parts
#   - .gs_output_dir()              The output folder of a part
#   - .gs_seed_example()            The example's pre-built files, put in place
#   - .gs_seed_closure(), .gs_seed_settings()   The settings the pre-built files stand for
#   - .gs_run_new(), .gs_value(), .gs_ancestors()        The state of one run
#   - .gs_unmet()                   The conditions of the table of parts
#   - .gs_part_folder(), .gs_unit_name(), .gs_expected_units(), .gs_unit_args()   Parts and units
#   - .gs_gene_groups()             The gene groups of the native coding comparison
#   - .gs_setting_keys(), .gs_setting_text()    The settings of a part, as text
#   - .gs_selected_files()          The files that are read from a folder
#   - .gs_unit_dependencies()       What a unit depends on
#   - .gs_locate()                  The analysis that holds a finished part
#   - .gs_folder_reusable()         Whether a whole folder is finished up the chain
#   - .gs_check_inputs()            The missing-input and example-data checks
#   - .gs_write_resolved()          config-resolved.yml of a step folder
#   - .gs_finish_early()            End a script whose unit has nothing to do
#
# in the GLOWpipeline monorepo.

#################### MAIN FUNCTIONS ####################

# The start of a stage script.
#
# step     "data_prep" | "training" | "single_variant" | "snv_set".
# part     The part of the step that the script is (for example "b_estimate"). It
#          decides the output folder.
# options  Named character vector of the script's own command-line options: the
#          name without the dashes, and "value" (the option takes a value) or "flag".
# load     TRUE: the script reads settings, so they are loaded into `envir`. FALSE:
#          the script needs its output folder only (the scan and summary scripts of
#          the SNV-set step read the snapshot that the preparation wrote).
# part_if  Named character vector: a flag of the script and the part that the
#          script is when the flag is set. marginal-scan.R is the part `scan`, and
#          with --combine the part `scan_combine`.
# records  TRUE: the unit of work is checked, locked and recorded. FALSE: the
#          settings and folders are given and nothing is
#          checked or written, for a tool that only looks at a configuration.
# envir    The environment of the script. Settings are assigned here.
# args     The command line.
# script   The file of the calling script, as .gs_running_script() finds it. It has to be
#          the file that base-config/parts.yml names for the part (and, for a scan, for
#          the region type), compared after links are resolved. NA makes no check, for a
#          call that is not a stage script. "" (the process does not say which script
#          it runs) is refused.
#
# Returns the context of the run, a list:
#   config_path  the configuration file
#   opts         the script's own options as a named list (a flag is TRUE or FALSE)
#   output_dir   where the part writes
#   run_name     the name of the run in logs, plot titles and run tables
#   sv_stem      the stem of the single-variant file names, "single_variant"
#   region_type  the region type of an SNV-set run, else NULL
#   cfg, layout  the resolved configuration and the layout
#   paths        (data_prep) the folders of the data-preparation parts
#   part, unit   the part and the unit of work that this process is
#   record       the state of the unit's completion record, for
#                glow_stage_end()
#   skipped      NULL, or why there is nothing to do: "complete" or
#                "reused". The script then ends at once with status 0.
#
# Assumption: the working directory is the GLOWanalyses directory, as for every
# stage script.
#
# The unit of work is checked before the script goes on
# (see the plan 16-config-records-plan.md). In this order: a folder in which every
# part of this run is finished in an analysis up the `extends` chain is reused, and
# the script ends.
# An input that no analysis holds stops the run and names the part to run first. An
# input inside the example data stops a run of another dataset. A unit whose record
# states the same dependencies is complete, and the script ends. A record with other
# dependencies stops the run unless --force is given.
glow_stage_begin <- function(step, part, options = character(0), load = TRUE, part_if = NULL,
                             records = TRUE, envir = parent.frame(), args = commandArgs(trailingOnly = TRUE),
                             script = .gs_running_script()) {
  steps <- c("data_prep", "training", "single_variant", "snv_set")
  if (!step %in% steps) stop("glow_stage_begin(): unknown step '", step, "'.")
  pa <- .gs_parse_args(args, options)

  # ---- the form of the configuration ----
  config_path <- pa$config
  if (is.na(config_path)) config_path <- "example.yml"
  if (!file.exists(config_path))
    stop(sprintf("Configuration not found: %s  (working directory: %s)", config_path, getwd()),
         call. = FALSE)
  # Release 0.3.6 took a config.R, which sourced R base files. One YAML file replaces it.
  if (grepl("\\.[Rr]$", config_path))
    stop("`", config_path, "` is a configuration in the form of release 0.3.6 (a config.R). GLOWanalyses reads one ",
         "YAML configuration: see MIGRATION.md in the GLOWanalyses directory for where each field of a config.R goes.",
         call. = FALSE)
  if (!grepl("\\.ya?ml$", config_path, ignore.case = TRUE))
    stop("The configuration has to be a YAML file (.yml): ", config_path, call. = FALSE)

  st <- list(step = step, part = part, config_path = config_path, opts = pa$opts,
             region_type = NULL, cfg = NULL, layout = NULL, paths = NULL, unit = NULL, record = NULL,
             skipped = NULL)

  # A value that release 0.3.6 took as an option is a setting here. It is refused before
  # anything is locked or written (see glow_stage_option()).
  for (o in setdiff(names(pa$opts), c("chr", "chunk", "arrayid", "spa", "combine")))
    if (!is.null(pa$opts[[o]]) && !identical(pa$opts[[o]], FALSE)) .gs_refuse_option(o)
  for (f in c("base-config/base_config.yml", "base-config/settings-registry.yml", "base-config/parts.yml"))
    if (!file.exists(f))
      stop("Run from the GLOWanalyses directory: ", f, " was not found in ", getwd(), ".", call. = FALSE)
  if (!requireNamespace("GLOWpipeline", quietly = TRUE) ||
      !exists("glow_unit_begin", envir = asNamespace("GLOWpipeline"), inherits = FALSE))
    stop("GLOWanalyses needs a GLOWpipeline version that has glow_layout() and glow_unit_begin().",
         call. = FALSE)

  cfg <- GLOWpipeline::read_glow_config(config_path, base = "base-config/base_config.yml",
                                        registry = "base-config/settings-registry.yml", root = ".")
  cfg <- GLOWpipeline::glow_absolute_paths(cfg)
  .gs_check_names(cfg)
  lay <- GLOWpipeline::glow_layout(cfg)
  if (step == "snv_set") st$region_type <- .gs_region_type(cfg, pa$region_type)
  else if (!is.null(pa$region_type)) stop("--region-type belongs to the SNV-set step.", call. = FALSE)

  # ---- the part and the unit of work that this process is ----
  parts <- GLOWpipeline::glow_read_parts("base-config/parts.yml", "base-config/settings-registry.yml")
  part_now <- part
  for (flag in names(part_if)) if (isTRUE(pa$opts[[flag]])) part_now <- part_if[[flag]]
  spec <- parts[[paste(step, part_now, sep = ".")]]
  if (is.null(spec)) stop("base-config/parts.yml has no part `", part_now, "` in the step `", step, "`.", call. = FALSE)
  # The script that runs has to be the script of this part. The four scan scripts share one
  # part, and each computes its own region type whatever the configuration says: without this
  # check, run-gene.R with a configuration of the coding run wrote gene tables into the coding
  # folder and completed the coding unit there. Nothing is seeded, locked or written before it.
  # The two are compared as files, with links resolved, and not by their names: a link or a
  # copy named run-coding.R that holds the gene script is not the coding script.
  expected_script <- sub("{region_type}", if (is.null(st$region_type)) "" else st$region_type, spec$script, fixed = TRUE)
  if (!is.na(script)) .gs_check_script(script, expected_script, step, part_now, st$region_type)
  .gs_check_control_source(cfg, spec)
  run <- .gs_run_new(cfg, lay, parts, st$region_type)

  st$layout     <- lay
  st$part       <- part_now
  st$output_dir <- .gs_output_dir(lay, step, part_now, st$region_type)
  st$run_name   <- paste(cfg$dataset, cfg$analysis, sep = "/")
  st$sv_stem    <- "single_variant"
  if (!records) {
    # a look at the configuration: the settings as the part would get them, nothing else
    st$cfg <- .gs_derive_inputs(run, spec)
    if (step == "data_prep") st$paths <- .gs_dataprep_paths(st$cfg, lay)
    if (load) .gs_assign_settings(st$cfg, step, st$region_type, envir)
    return(st)
  }
  if (identical(cfg$dataset, "example")) .gs_seed_example(run)
  unit    <- .gs_unit_name(spec, pa$opts, st$region_type)
  folder  <- .gs_part_folder(lay, spec, st$region_type)
  what    <- paste0(step, "/", if (is.null(unit)) part_now else paste(part_now, unit))
  st$unit <- unit
  # The units of a step folder: the record of the folder holds the derived inputs of every
  # part that writes there, so it does not depend on which part wrote it last.
  same_dir <- Filter(function(p) identical(p$step, step) &&
                       identical(.gs_output_dir(lay, step, p$part, st$region_type), st$output_dir), parts)

  if (identical(cfg$run$records, FALSE)) {
    # A workflow engine runs this unit in a container of its own and keeps track of what is
    # finished (the setting run.records). The unit gets its settings and its folders and
    # writes the resolved configuration. It keeps no record, takes no lock and checks no
    # input by a record: in a new container there is nothing to compare with.
    st$cfg <- .gs_derive_inputs(run, spec)
    if (step == "data_prep") st$paths <- .gs_dataprep_paths(st$cfg, lay)
    .gs_write_resolved(run, .gs_derive_inputs(run, unname(same_dir)), st$output_dir, what)
    if (load) .gs_assign_settings(st$cfg, step, st$region_type, envir)
    return(st)
  }

  # ---- the checks before any work, then the lock and the overwrite guard ----
  deps <- .gs_unit_dependencies(run, spec$id, unit, note = TRUE)
  if (!pa$force && spec$level != "dataset" && .gs_folder_untouched(folder)) {
    from <- .gs_folder_reusable(run, spec)
    if (!is.null(from)) {
      st$cfg <- cfg
      return(.gs_finish_early(st, "reused", sprintf(
        "%s: reused from the analysis `%s`, whose results in %s were made with the same dependencies. Nothing to do.",
        what, from$analysis, from$folder)))
    }
  }
  .gs_check_inputs(run, spec, unit, deps)
  handle <- GLOWpipeline::glow_unit_begin(folder, spec$part, unit, deps, force = pa$force,
                                          level = if (spec$level == "dataset") "dataset" else "analysis")
  if (identical(handle$status, "complete")) {
    st$cfg <- cfg
    return(.gs_finish_early(st, "complete", sprintf("%s: complete (record %s). Nothing to do.",
                                                    what, handle$record_file)))
  }
  if (length(handle$stale))
    cat(what, ": runs again, because a part it read was produced again since it ran:\n",
        paste0("  - ", handle$stale, collapse = "\n"), "\n", sep = "")
  # From here the unit runs and holds its lock. The lock is released by
  # glow_stage_end(), or when R exits without it (an error in the script).
  state <- new.env()
  state$handle <- handle; state$done <- FALSE
  state$spec <- spec; state$view <- NULL
  state$run_new <- function() .gs_run_new(cfg, lay, parts, st$region_type)
  state$deps_now <- function() .gs_unit_dependencies(.gs_run_new(cfg, lay, parts, st$region_type), spec$id, unit)
  reg.finalizer(state, function(e) {
    if (!is.null(e$view)) unlink(e$view, recursive = TRUE, expand = FALSE)            # links only, see glow_stage_input_dir()
    if (!e$done) GLOWpipeline::glow_unit_release(e$handle)
  }, onexit = TRUE)
  st$record <- state

  # ---- the settings as the part uses them, and the record of them ----
  # (`cfg` stays the configuration as the user wrote it: the dependencies are read from it.)
  cfg_used <- .gs_derive_inputs(run, spec)
  st$cfg <- cfg_used
  if (step == "data_prep") st$paths <- .gs_dataprep_paths(cfg_used, lay)
  .gs_write_resolved(run, .gs_derive_inputs(run, unname(same_dir), note = TRUE), st$output_dir, what)
  if (load) .gs_assign_settings(cfg_used, step, st$region_type, envir)
  st
}

# The state of a unit, without running it.
#
# The `glow` command asks this for every unit of a step: to show what a call would
# do (--dry-run), and to decide which units to submit as jobs. Nothing is locked or
# written.
#
# cfg, lay     The configuration (paths absolute) and its layout.
# parts        The table of parts.
# spec         The part.
# region_type  The region type of an SNV-set part, else NULL.
# unit         The unit, or NULL for a part that is not split.
# force        Whether the call is forced.
# run          The state of a run (.gs_run_new()), to share what was looked up between
#              the units of one call. NULL: a new one.
#
# Returns list(status, detail). `status` is one of
#   complete  the record states the same dependencies and its outputs exist
#   reused    the folder is finished in an analysis up the `extends` chain
#   ready     the unit can run now
#   waiting   a part it reads is not finished (`detail` names the parts)
#   refused   the folder holds results made with other settings; --force replaces them
glow_unit_status <- function(cfg, lay, parts, spec, region_type, unit, force = FALSE, run = NULL) {
  if (is.null(run)) run <- .gs_run_new(cfg, lay, parts, region_type)
  run$region_type <- region_type
  folder <- .gs_part_folder(lay, spec, region_type)
  deps <- .gs_unit_dependencies(run, spec$id, unit)
  if (!force && spec$level != "dataset" && .gs_folder_untouched(folder)) {
    from <- .gs_folder_reusable(run, spec)
    if (!is.null(from)) return(list(status = "reused", detail = paste("from the analysis", from$analysis)))
  }
  missing <- names(deps$inputs)[deps$inputs == "missing"]
  state <- GLOWpipeline::glow_unit_state(folder, spec$part, unit, deps)
  if (length(missing))
    return(list(status = "waiting", detail = paste("for", paste(sub("^[a-z_]+\\.", "", missing), collapse = ", ")),
                missing = unique(sub("^([a-z_]+\\.[a-z_]+)\\..*$", "\\1", missing))))
  if (force) return(list(status = "ready", detail = "forced"))
  if (state$complete) return(list(status = "complete", detail = ""))
  if (!is.null(state$record) && length(state$differences) && !isTRUE(state$stale))
    return(list(status = "refused", detail = paste(state$differences, collapse = "; ")))
  list(status = "ready", detail = if (is.null(state$record)) "" else state$reason)
}

# The end of a stage script: the completion record of its unit.
#
# st       The context that glow_stage_begin() returned.
# outputs  The files and folders the unit produced. Each has to exist.
#
# The record is written and the lock of the unit is released (nothing happens for a
# run without records). The inputs and the outside files are looked up
# once more first: when one of them differs from the begin of the unit, no record
# is written and the script stops, because the results belong to an input that is
# gone. A folder among the outputs stands for every file in it.
#
# Returns the record, invisibly.
glow_stage_end <- function(st, outputs = character(0)) {
  if (is.null(st$record)) return(invisible(NULL))
  state <- st$record
  state$done <- TRUE     # glow_unit_complete() releases the lock, also when it stops
  if (!is.null(state$view)) { unlink(state$view, recursive = TRUE, expand = FALSE); state$view <- NULL }
  now <- state$deps_now()
  record <- GLOWpipeline::glow_unit_complete(state$handle, outputs = outputs, inputs_now = now$inputs,
                                             files_now = now$files, producer = .gs_producer())
  invisible(record)
}

# Which process produces the unit, for its record: the job of the scheduler, when the
# script runs in one. `glow run --wait` reads it to tell whether a record is the one that
# the job it followed wrote. NULL for a run outside a job.
.gs_producer <- function() {
  job <- Sys.getenv("SLURM_JOB_ID", "")
  if (nzchar(job)) c(slurm_job = job) else NULL
}

# The finished results of this run that a gathering part reads.
#
# The aggregation gathers the tables of the scan units. The results folder can hold
# more than the tables of this run: the table of a unit that is not finished, or of a
# unit that belonged to an earlier configuration of the same folder (a chromosome that
# was taken out). A part must read the results that its completion record covers, and
# no other. So this function builds a view: a folder of links
# to the files that the completion records of the finished units list, for every part
# of the same folder that the running part reads, under the same relative names. The
# part reads from the view and writes to its real folder. With `allow_incomplete`, a
# unit that is not finished is simply not in the view.
#
# st  The context of the run.
#
# Returns the folder to read from: the view for a run that keeps records, and
# st$output_dir otherwise (a task of a workflow engine,
# whose results folder holds what the engine put there). glow_stage_end() removes the
# view. It holds links only.
glow_stage_input_dir <- function(st) {
  state <- st$record
  if (is.null(state)) return(st$output_dir)
  run <- state$run_new(); spec <- state$spec
  folder <- .gs_part_folder(run$layout, spec, run$region_type)
  files <- character(0)
  for (inp in spec$inputs) {
    other <- run$parts[[inp$part]]
    if (!identical(other$folder, spec$folder) || !identical(other$level, spec$level)) next
    if (nzchar(.gs_unmet(inp$when, run$cfg, run$region_type))) next
    for (u in .gs_expected_units(run, other, folder)) {
      unit_state <- GLOWpipeline::glow_unit_state(folder, other$part, u, .gs_unit_dependencies(run, other$id, u))
      if (unit_state$complete) files <- c(files, unit_state$record$outputs)
    }
  }
  files <- unique(files[!grepl("^/", files)])
  view <- file.path(folder, paste0(".inputs-", spec$part))
  unlink(view, recursive = TRUE, expand = FALSE)     # of a run that did not reach its end
  dir.create(view, recursive = TRUE)
  state$view <- view
  for (rel in files) {
    dir.create(dirname(file.path(view, rel)), recursive = TRUE, showWarnings = FALSE)
    file.symlink(file.path(folder, rel), file.path(view, rel))
  }
  view
}

# The files that a part wrote into a folder it shares with another part.
#
# The aggregation and the plots of an SNV-set run write into one folder. The plots
# record the files they wrote themselves: glow_stage_watch() looks at the folder before
# the part writes, and glow_stage_written() gives the files that are new or changed
# since, hidden files included. This holds while no other part writes into the folder
# at the same time. When one does (an aggregation that is forced while the plots run),
# the plots do not complete: their input was replaced, and the completion check says
# so. The aggregation lists its own tables by name and does not use the watcher.
#
# dir     The folder.
# before  What glow_stage_watch() returned.
glow_stage_watch <- function(dir) {
  files <- list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  info <- file.info(files)
  stats::setNames(sprintf("%.0f@%.6f", info$size, as.numeric(info$mtime)), files)
}
glow_stage_written <- function(dir, before) {
  now <- glow_stage_watch(dir)
  names(now)[!(names(now) %in% names(before)) | now != before[names(now)]]
}

# An empty folder for a part that owns one.
#
# A part that owns a folder (the PI models) is recorded with every file in it, and
# the next part reads every file in it. So the folder is emptied before the part
# writes: a file that a run left there before it stopped must not pass as a result
# of this run. Only for a folder inside the
# part's own output folder. A link in the folder is removed as a link: what it
# points at is not touched, and a linked folder is not entered.
#
# st   The context of the run.
# dir  The folder.
glow_stage_clear_dir <- function(st, dir) {
  if (!dir.exists(dir)) return(invisible(NULL))
  inside <- startsWith(normalizePath(dir), paste0(normalizePath(st$output_dir), "/"))
  if (!inside) stop("glow_stage_clear_dir(): ", dir, " is not inside the output folder of the part.")
  if (basename(dir) %in% .gs_reserved_names) stop("glow_stage_clear_dir(): ", dir, " is a folder of the workflow itself.")
  empty <- function(d) for (p in list.files(d, full.names = TRUE, all.files = TRUE, no.. = TRUE)) {
    if (nzchar(Sys.readlink(p)) || !dir.exists(p)) unlink(p, expand = FALSE) else empty(p)
  }
  empty(dir)
  invisible(NULL)
}

# The files of one unit of an SNV-set scan.
#
# GLOWpipeline::run_scan_unit() writes the unit's tables into the results folder
# (the table, and beside it the errors and the summary by category) and, when the
# output policy asks, one file per unit in evidence/, staar_detail/ and
# result_objects/. A script calls this function twice: with clear = TRUE before the
# scan, which removes what an earlier run of the unit left, and at the end, for the
# completion record. So the record lists the files of this run and all of them.
#
# st           The context of the run.
# results_dir  The results folder of the scan.
# stem         The name of the unit's table without its extension (glow_chr21).
# unit_id      The unit as run_scan_unit() names it (chr21, chunk_0003).
# clear        TRUE: delete the files and return nothing.
glow_scan_unit_files <- function(st, results_dir, stem, unit_id, clear = FALSE) {
  kinds <- c("evidence", "staar_detail", "result_objects")
  tables <- list.files(results_dir, full.names = TRUE)
  tables <- tables[startsWith(basename(tables), paste0(stem, ".")) | startsWith(basename(tables), paste0(stem, "_"))]
  side <- file.path(st$output_dir, kinds, sprintf("%s_%s.rds", kinds, unit_id))
  files <- c(tables, side[file.exists(side)])
  if (!clear) return(files)
  unlink(files, expand = FALSE)
  invisible(NULL)
}

# Links from an analysis to the data of its dataset.
#
# The GDS files, the annotated GDS files and the principal components belong to the
# dataset and are written once under <dataset>/data/. The data_prep folder of an
# analysis holds the phenotype bundle and one relative link per data folder that
# exists, so that a reader of the analysis finds the data it used.
#
# lay            The layout of the run (GLOWpipeline::glow_layout()).
# data_prep_dir  The data_prep folder that gets the links. Default: the one of the
#                run's own analysis.
#
# Returns the paths of the links, invisibly. Deleting the analysis folder deletes
# the links only. `rm -rf <link>/`, with the trailing slash, deletes the data.
glow_link_dataset_assets <- function(lay, data_prep_dir = lay$steps$data_prep) {
  is_link <- function(p) { target <- Sys.readlink(p); !is.na(target) && nzchar(target) }
  dir.create(data_prep_dir, recursive = TRUE, showWarnings = FALSE)
  links <- character(0)
  for (name in c("gds", "gds_favor", "pcs")) {
    link <- file.path(data_prep_dir, name)
    # two levels up from <analysis>/data_prep/ is the dataset folder
    if (dir.exists(file.path(lay$data$dir, name)) && !file.exists(link) && !is_link(link))
      file.symlink(file.path("..", "..", "data", name), link)
    if (is_link(link)) links <- c(links, link)
  }
  invisible(links)
}

#################### INTERNAL HELPERS ####################

# The command line of a stage script.
#
# args     Character vector, as commandArgs(trailingOnly = TRUE) gives it.
# options  The script's own options (see glow_stage_begin()).
#
# Returns list(config, region_type, force, opts). `config` is NA when --config is
# absent. An argument that is neither --config, --region-type, --force nor a
# declared option stops the script, as the scripts of release 0.3.6 do.
.gs_parse_args <- function(args, options) {
  config <- NA_character_; region_type <- NULL; force <- FALSE
  opts <- lapply(options, function(kind) if (kind == "flag") FALSE else NULL)
  names(opts) <- names(options)
  i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    name <- sub("^--", "", a)
    has_value <- i < length(args)
    if (a == "--config" && has_value) { config <- args[i + 1L]; i <- i + 2L }
    else if (a == "--region-type" && has_value) { region_type <- args[i + 1L]; i <- i + 2L }
    else if (a == "--force") { force <- TRUE; i <- i + 1L }
    else if (startsWith(a, "--") && name %in% names(options) && options[[name]] == "flag") {
      opts[[name]] <- TRUE; i <- i + 1L
    } else if (startsWith(a, "--") && name %in% names(options) && has_value) {
      opts[name] <- list(args[i + 1L]); i <- i + 2L
    } else stop("Unknown argument: ", a, call. = FALSE)
  }
  list(config = config, region_type = region_type, force = force, opts = opts)
}

# The one region type of an SNV-set run.
#
# cfg       The resolved configuration. `snv_set.region_type` lists the types to run.
# selected  The value of --region-type, or NULL.
#
# Returns the type. A list of one type needs no selector. A selector has to name a
# listed type.
# A part that draws the PI controls needs a folder of annotated per-chromosome files.
#
# By default that is the folder the annotation of this dataset writes. A study that
# supplies its annotated GDS files (shared.analysis_gds_dir) leaves the annotation out
# of the run, so nothing writes the default folder, and the part would wait for a part
# that cannot run. The folder is then the user's to name. It is not taken from the
# supplied GDS files by itself: the PI training finds control files by the names
# chr<N>...gds or chr<N>...csv, which supplied files need not have.
#
# cfg   The configuration as the user wrote it (paths absolute).
# spec  The part.
.gs_check_control_source <- function(cfg, spec) {
  reads <- any(vapply(spec$inputs, function(inp) identical(inp$setting, "training.pi_control_source"), logical(1)))
  if (!reads || !is.null(cfg$training$pi_control_source) || is.null(cfg$shared$analysis_gds_dir)) return(invisible(NULL))
  stop("`", spec$step, "/", spec$part, "` draws the control variants of the PI models from annotated per-chromosome files. ",
       "This configuration supplies the annotated GDS files (`shared.analysis_gds_dir`), so the data preparation writes ",
       "no such files and `training.pi_control_source` has to name their folder. The folder of the supplied GDS files ",
       "serves when its files are named chr<N>...gds: set `training.pi_control_source` to it and ",
       "`training.pi_control_format` to `gds`.", call. = FALSE)
}

# The file of the R script that runs now, or "" when the process does not say.
#
# A script that is being read by source() is the innermost such file. Otherwise it is the
# file that R was started with (`Rscript <file>`, `R -f <file>`). "" is not NA: NA is what
# a caller passes to glow_stage_begin() to make no check (a test, a look at a
# configuration), and "" is a process whose script cannot be told, which is refused.
.gs_running_script <- function() {
  for (i in rev(seq_len(sys.nframe()))) {
    if (!identical(sys.function(i), base::source)) next
    of <- get0("ofile", envir = sys.frame(i), inherits = FALSE)
    if (is.character(of) && length(of) == 1L && nzchar(of)) return(of)
  }
  args <- commandArgs(FALSE)
  file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (!length(file)) {
    f <- which(args == "-f")
    if (length(f) && f[1L] < length(args)) file <- args[f[1L] + 1L]
  }
  if (length(file) == 1L && nzchar(file)) file else ""
}

# Stop unless `script` is the file of the part's script.
#
# script    The file that runs, as .gs_running_script() gives it ("" when unknown).
# expected  The script of the part in base-config/parts.yml, relative to the GLOWanalyses
#           directory, which is the working directory.
# The comparison is between the two files after links are resolved, so a link to the
# right script passes under any name, and a link or a copy that carries the right name
# and another script does not.
.gs_check_script <- function(script, expected, step, part, region_type) {
  what <- paste0("the part `", step, "/", part, "`",
                 if (!is.null(region_type)) paste0(" for the region type `", region_type, "`") else "")
  hint <- paste0(" `glow run ", step, " --part ", part, "` starts the script that belongs to the configuration.")
  if (!nzchar(script))
    stop("This R process does not say which script it runs: it was not started as `Rscript <script>` and the ",
         "script is not read by `source()`. So the script cannot be checked against ", what, ", whose script is `",
         expected, "`.", hint, call. = FALSE)
  same <- file.exists(script) && file.exists(expected) &&
    identical(normalizePath(script, mustWork = TRUE), normalizePath(expected, mustWork = TRUE))
  if (!same)
    stop("`", script, "` is not the script of ", what, ": that is `", expected, "`.", hint, call. = FALSE)
  invisible(NULL)
}

.gs_region_type <- function(cfg, selected) {
  listed <- cfg$snv_set$region_type
  if (is.null(selected)) {
    if (length(listed) == 1L) return(listed)
    stop("The configuration lists ", length(listed), " region types (", paste(listed, collapse = ", "),
         "). This script runs one: pass --region-type <type>, as `glow run snv_set` does.", call. = FALSE)
  }
  if (!selected %in% listed)
    stop("--region-type ", selected, " is not among the region types of the configuration (",
         paste(listed, collapse = ", "), "). Add it to `snv_set.region_type`.", call. = FALSE)
  selected
}

# Where each derived input comes from.
#
# A step finds the outputs of the parts before it through the layout, so a user
# writes no path between steps. Each rule below belongs to one setting whose
# default is null:
#   part       the part that produces the file or folder (step.part), or NULL for a
#              value that is not a path
#   producers  the parts for which the setting is where they write, not what they
#              read. Default: `part`.
#   when       a condition on the configuration without which the setting has no
#              place to stand in (a field of a map that is null)
#   path       the place in a layout, as function(lay, cfg, run)
#
# Whether a part reads an input in a run is said in the table of parts, by the
# `when` of the input, and nowhere else.
#
# Returns the rules, a list named by setting (section.name, or section.name.field
# for a field of a map).
.gs_derived_rules <- function() {
  list(
    "data_prep.favor_input_gds_dir" = list(part = "data_prep.plink_to_gds",
      path = function(lay, cfg, run) lay$data$gds),
    "data_prep.pcs.path" = list(part = "data_prep.compute_pcs",
      when = function(cfg) !is.null(cfg$data_prep$pcs),
      path = function(lay, cfg, run) file.path(lay$data$pcs, "pcs.rds")),
    "shared.analysis_gds_dir" = list(part = "data_prep.annotate_favor",
      path = function(lay, cfg, run) lay$data$gds_favor_gds),
    "shared.analysis_gds_pattern" = list(part = NULL,
      path = function(lay, cfg, run) cfg$data_prep$favor_gds_pattern),
    "shared.pheno_path" = list(part = "data_prep.assemble_pheno",
      path = function(lay, cfg, run) file.path(lay$steps$data_prep, cfg$data_prep$pheno_out_name)),
    "training.b_training_rds" = list(part = "training.b_prepare",
      path = function(lay, cfg, run) file.path(lay$steps$training_b, "b_training_data.rds")),
    "training.pi_case_csv" = list(part = "training.pi_prepare",
      path = function(lay, cfg, run) file.path(lay$steps$training_pi, "pi_case_data.csv")),
    "training.pi_case_annotated_csv" = list(part = "training.pi_annotate",
      path = function(lay, cfg, run) file.path(lay$steps$training_pi, "pi_case_data_annotated.csv")),
    "training.pi_models_dir" = list(part = "training.pi_train",
      path = function(lay, cfg, run) file.path(lay$steps$training_pi, cfg$training$pi_models_subdir)),
    "training.pi_control_source" = list(part = "data_prep.annotate_favor",
      path = function(lay, cfg, run) .gs_control_source(run, lay)),
    "single_variant.ld_scores_path" = list(part = "single_variant.ld_scores_combine",
      producers = c("single_variant.ld_scores", "single_variant.ld_scores_combine"),
      path = function(lay, cfg, run) file.path(lay$steps$single_variant, "ld_scores", "ld_scores.rds")),
    "snv_set.pi_model_dir" = list(part = "training.pi_train",
      path = function(lay, cfg, run) file.path(lay$steps$training_pi, cfg$training$pi_models_subdir)),
    "snv_set.b_model_path" = list(part = "training.b_estimate",
      path = function(lay, cfg, run) file.path(lay$steps$training_b, "b_model.rds")),
    "snv_set.calibration.report" = list(part = "single_variant.calibrate",
      path = function(lay, cfg, run) file.path(lay$steps$single_variant, "results", "single_variant_calibration.tsv")),
    "snv_set.calibration.ld_scores" = list(part = "single_variant.ld_scores_combine",
      path = function(lay, cfg, run) file.path(lay$steps$single_variant, "ld_scores", "ld_scores.rds")))
}

# Inputs left null, filled from the layout.
#
# An input that the configuration sets is used as written. An input left null is
# filled here for the part that runs, and only when the part reads it in this run
# (the `when` of the input in the table of parts). A setting that the part writes
# to is placed in the run's own analysis. A setting that the part reads is placed
# in the analysis where the producing part was found finished: the run's own, or
# one up the `extends` chain. Every other derived input stays null, because this
# part does not use it. A choice that the run makes for the part (the control
# format that "auto" becomes) is filled in the same way, so that the script gets
# the value that the record states.
#
# run    The state of the run (.gs_run_new()).
# specs  The part that runs, from the table of parts, or a list of parts: the parts
#        that write into one step folder, for the record of that folder.
# note   Whether an input found up the chain is noted as reused.
#
# Returns the configuration with the derived inputs of these parts filled in.
.gs_derive_inputs <- function(run, specs, note = FALSE) {
  if (!is.null(specs$id)) specs <- list(specs)
  cfg <- run$cfg; rules <- .gs_derived_rules()
  ids <- vapply(specs, `[[`, character(1), "id")
  # the parts that read the setting in this run. The file-name pattern goes with the
  # folder it describes.
  readers_of <- function(key) Filter(function(spec) any(vapply(spec$inputs, function(inp)
    (identical(inp$setting, key) ||
       (key == "shared.analysis_gds_pattern" && identical(inp$setting, "shared.analysis_gds_dir"))) &&
      !nzchar(.gs_unmet(inp$when, run$cfg, run$region_type)), logical(1))), specs)
  for (key in names(rules)) {
    rule <- rules[[key]]
    producers <- if (is.null(rule$producers)) rule$part else rule$producers
    writes <- any(ids %in% producers)
    readers <- readers_of(key)
    if (!writes && !length(readers)) next
    if (!is.null(.gs_value(cfg, key))) next                       # set by the user
    if (!is.null(rule$when) && !isTRUE(rule$when(run$cfg))) next   # the setting has no place
    lay <- run$layout
    if (!writes && !is.null(rule$part)) {
      found <- .gs_locate(run, rule$part, from = readers[[1L]], note = note)
      if (isTRUE(found$found)) lay <- found$layout
    }
    cfg <- .gs_set_value(cfg, key, rule$path(lay, run$cfg, run))
  }
  for (key in "training.pi_control_format")
    if (any(vapply(specs, function(spec) key %in% spec$settings, logical(1))))
      cfg <- .gs_set_value(cfg, key, .gs_effective_value(run, key))
  cfg
}

# A setting's value as the run uses it.
#
# Two settings leave a choice open that the run makes. A null file-name pattern of
# the analysis GDS files is the pattern of the annotation step. The control format
# "auto" of the PI training, with derived control files, is "csv" or "gds"
# (.gs_control_format()). The completion record states this value, and the script
# gets it, so a changed choice is a changed dependency.
.gs_effective_value <- function(run, key) {
  cfg <- run$cfg
  value <- .gs_value(cfg, key)
  if (key == "shared.analysis_gds_pattern" && is.null(value)) return(cfg$data_prep$favor_gds_pattern)
  if (key == "training.pi_control_format" && is.null(cfg$training$pi_control_source)) return(.gs_control_format(run))
  value
}

# The format of the derived control files of the PI training: "csv" or "gds".
#
# `pi_control_format` names it, and "agds" is "gds". "auto" is "csv" when the
# annotation step wrote a CSV file for every chromosome of `pi_chromosomes`: the
# file exists and the completion record of that chromosome's annotation lists it.
# Otherwise it is "gds". A CSV file that was put in the folder by hand is not what
# "auto" reads.
.gs_control_format <- function(run) {
  cfg <- run$cfg; lay <- run$layout
  fmt <- tolower(cfg$training$pi_control_format)
  if (fmt == "csv") return("csv")
  if (fmt %in% c("gds", "agds")) return("gds")
  chroms <- cfg$training$pi_chromosomes
  if (!length(chroms)) return("gds")
  folder <- lay$data$gds_favor                       # the records of the annotation are kept here
  for (k in chroms) {
    csv <- file.path(lay$data$gds_favor_csv, gsub("{chr}", k, cfg$data_prep$favor_csv_pattern, fixed = TRUE))
    rec <- GLOWpipeline::glow_unit_record(folder, "annotate_favor", paste0("chr", k))
    listed <- !is.null(rec) && substring(csv, nchar(folder) + 2L) %in% rec$outputs   # as the record names it
    if (!listed || !file.exists(csv)) return("gds")
  }
  "csv"
}

# The folder of the derived control files, in a layout: its csv/ or gds/ child.
.gs_control_source <- function(run, lay)
  if (.gs_control_format(run) == "csv") lay$data$gds_favor_csv else lay$data$gds_favor_gds

# Annotation tests as GLOWr objects.
#
# A configuration writes a test of an annotation as a map with `op` and, for most
# operators, `value`. GLOWr::variant_filter() and the coding categories take the
# object that GLOWr::annotation_predicate() builds. This function walks a value
# (the clauses of a filter, or custom categories) and converts every such map.
#
# x  A list of clauses, a list of categories, or any part of one.
#
# Returns x with every test converted. Everything else is returned as it is.
.gs_as_predicates <- function(x) {
  ops <- c("in", "not_in", "gt", "ge", "lt", "le", "nonempty", "empty")
  if (!is.list(x) || inherits(x, "glow_annotation_predicate")) return(x)
  is_test <- !is.null(names(x)) && all(names(x) %in% c("op", "value")) && "op" %in% names(x) &&
    is.character(x$op) && length(x$op) == 1L && x$op %in% ops
  if (is_test) return(GLOWr::annotation_predicate(x$op, x$value))
  for (i in seq_along(x)) x[i] <- list(.gs_as_predicates(x[[i]]))
  x
}

# The settings under the scripts' variable names.
#
# Each registry entry records the base file and the name its setting had in release
# 0.3.6 (`was`). For the step, every setting whose `was` names the step's base file
# is assigned under that old name, as the R object the scripts expect.
#
# cfg          The resolved configuration, derived inputs filled in.
# step         The step of the script.
# region_type  The one region type of an SNV-set run, else NULL.
# envir        The environment of the script.
#
# Returns the assigned names, invisibly.
.gs_assign_settings <- function(cfg, step, region_type, envir) {
  base_file <- c(data_prep = "data_prep_base.R", training = "training_base.R",
                 single_variant = "single_variant_base.R", snv_set = "snv_set_base.R")[[step]]
  assigned <- character(0)
  for (e in attr(cfg, "glow_registry")$entries) {
    if (is.null(e$was) || !base_file %in% unlist(e$was$file)) next
    old_name <- if (is.null(e$was$name)) e$name else e$was$name
    value <- if (e$section == "top") cfg[[e$name]] else cfg[[e$section]][[e$name]]
    # an empty map is the plain empty list that the scripts expect (list())
    if (is.list(value) && length(value) == 0L) value <- list()
    assign(old_name, value, envir = envir)
    assigned <- c(assigned, old_name)
  }
  # GLOWpipeline::prepare_scan_run() links the outputs into a lab's shared storage
  # when this field is absent. NULL means a plain folder, which the layout places.
  assign("symlinked_shared_root", NULL, envir = envir)
  if (step == "snv_set") {
    assign("region_type", region_type, envir = envir)
    fs <- get("filter_spec", envir = envir)
    fs["annotation_clauses"] <- list(.gs_as_predicates(fs$annotation_clauses))
    assign("filter_spec", fs, envir = envir)
    assign("categories", .gs_as_predicates(get("categories", envir = envir)), envir = envir)
  }
  invisible(c(assigned, "symlinked_shared_root"))
}

# The folders of the data-preparation parts, under the names that
# resolve_cohort_paths() of 00-data-prep/_dataprep_lib.R uses. The three `read_`
# elements are what a part reads, which is not always what another part of the run
# wrote: the folder of GDS files and the folder and file-name pattern of the
# annotated GDS files are the user's when the configuration names them
# (data_prep.favor_input_gds_dir, shared.analysis_gds_dir and
# shared.analysis_gds_pattern), and the layout's otherwise. They are NULL for a part
# that does not read them.
.gs_dataprep_paths <- function(cfg, lay) {
  list(plink_dir       = cfg$data$plink_dir,
       gds_dir         = lay$data$gds,
       favor_dir       = lay$data$gds_favor,
       favor_input_dir = cfg$data_prep$favor_input_gds_dir,
       favor_gds_dir   = lay$data$gds_favor_gds,
       favor_csv_dir   = lay$data$gds_favor_csv,
       pcs_dir         = lay$data$pcs,
       pheno_dir       = lay$steps$data_prep,
       read_gds_dir           = cfg$data_prep$favor_input_gds_dir,
       read_favor_gds_dir     = cfg$shared$analysis_gds_dir,
       read_favor_gds_pattern = cfg$shared$analysis_gds_pattern)
}

# The output folder of a part. The data-preparation parts write to several
# folders and take them from st$paths, so their output_dir is the analysis's
# data_prep folder.
.gs_output_dir <- function(lay, step, part, region_type) {
  switch(step,
    data_prep      = lay$steps$data_prep,
    training       = if (startsWith(part, "b_")) lay$steps$training_b else lay$steps$training_pi,
    single_variant = lay$steps$single_variant,
    snv_set        = lay$steps$snv_set[[region_type]])
}

# The example's pre-built files, put in place, each unit with a seed record.
#
# The shipped example runs its single-variant and SNV-set steps on files that are
# built ahead: the annotated GDS files, the phenotype bundle, the B model and the PI
# models. The layout derives where a run reads them, so they are supplied there from
# data-example/, following data-example/seed-manifest.yml. The dataset folder gets
# the marker file EXAMPLE-DATA.
#
# Each seeded unit gets a seed record with the settings the files were built with.
# These are fixed in data-example/seed-settings.yml, which is written when the files
# are built, so a later edit of example.yml does not change what the files stand
# for. They cover the part itself and every part before it, with the outside files
# by their paths: a run that changes a setting anywhere up the chain is told to
# produce its own files. GLOWpipeline::glow_seed_unit() does the work: it holds the
# lock of the unit, never writes over the results of a run, and accepts a file
# that is in place only when it holds the supplied bytes.
#
# run  The state of the run (.gs_run_new()).
#
# Returns the status of each seeded unit, invisibly.
.gs_seed_example <- function(run) {
  manifest_file <- "data-example/seed-manifest.yml"; settings_file <- "data-example/seed-settings.yml"
  if (!file.exists(manifest_file)) return(invisible(character(0)))
  if (!file.exists(settings_file))
    stop("The example's pre-built files have no ", settings_file, ". Write it with data-example/write-seed-settings.R.",
         call. = FALSE)
  manifest <- yaml::read_yaml(manifest_file)
  frozen <- yaml::read_yaml(settings_file)
  problems <- .gs_seed_problems(manifest, frozen)
  if (length(problems))
    stop("The example's pre-built files under ", manifest$source, " are not the set that ", settings_file, " describes:\n",
         paste0("  - ", problems, collapse = "\n"),
         "\nRestore data-example/ as the repository has it (git checkout -- data-example), or rebuild the files with ",
         "data-example/build-intermediates.sh.", call. = FALSE)
  # the analysis that holds the seeded files of the analysis level, under this output root
  seed_cfg <- run$cfg; seed_cfg$analysis <- manifest$analysis
  seed_lay <- GLOWpipeline::glow_layout(seed_cfg)
  dataset_dir <- run$layout$dataset_dir
  units <- list()
  for (item in manifest$files) {
    from <- file.path(manifest$source, item$from)
    to   <- file.path(dataset_dir, item$to)
    if (!file.exists(from))
      stop("The example's seed manifest names a file that is missing: ", from, call. = FALSE)
    if (dir.exists(from)) {
      rel <- list.files(from, recursive = TRUE)
      from <- file.path(from, rel); to <- file.path(to, rel)
    }
    key <- paste(item$part, if (is.null(item$unit)) "" else item$unit)
    units[[key]] <- list(part = item$part, unit = item$unit, from = c(units[[key]]$from, from), to = c(units[[key]]$to, to))
  }
  status <- character(0)
  for (key in names(units)) {
    u <- units[[key]]
    spec <- run$parts[[u$part]]
    if (is.null(spec)) stop("The example's seed manifest names a part that base-config/parts.yml does not hold: ", u$part)
    built <- frozen[[u$part]]
    if (is.null(built)) stop(settings_file, " has no entry for the part ", u$part, ".", call. = FALSE)
    # the settings, and the content of the outside files, that the files were built with
    settings <- c(unlist(built$settings), unlist(built$inputs))
    status[[key]] <- GLOWpipeline::glow_seed_unit(.gs_part_folder(seed_lay, spec, NULL), spec$part, u$unit,
                                                  from = u$from, to = u$to, settings = settings)
  }
  marker <- file.path(dataset_dir, "EXAMPLE-DATA")
  if (!file.exists(marker) && dir.exists(dataset_dir))
    writeLines(c("This dataset folder belongs to the shipped GLOWanalyses example.",
                 "Its inputs come from data-example/, a small synthetic cohort.",
                 "Models and tables in this folder are valid for the example only."), marker)
  if (any(status == "seeded")) {
    cat(sprintf("Example: supplied the pre-built files of %d unit(s) from %s into %s\n",
                sum(status == "seeded"), manifest$source, dataset_dir))
    glow_link_dataset_assets(run$layout, seed_lay$steps$data_prep)
  }
  invisible(status)
}

# The marker of a collection of the pre-built files that has begun and not ended
# (data-example/write-seed-settings.R --collect).
GLOW_SEED_MARKER <- "data-example/seed-collection.incomplete"

# The pre-built files that are in place, by their paths under the folder that the
# manifest names. A folder of the manifest stands for every file in it.
.gs_seed_files <- function(manifest) {
  unlist(lapply(manifest$files, function(item) {
    from <- file.path(manifest$source, item$from)
    if (dir.exists(from)) file.path(item$from, list.files(from, recursive = TRUE, all.files = TRUE, no.. = TRUE)) else item$from
  }))
}

# What keeps the pre-built files from being used.
#
# The files and data-example/seed-settings.yml are written by one collection. The
# settings file states the checksum of every file (`artifacts`), so files of two
# collections, a file that was replaced and a file that is missing are noticed. A
# collection that was interrupted leaves GLOW_SEED_MARKER.
#
# manifest  The seed manifest, as read.
# frozen    The settings file, as read.
#
# Returns the problems, one string each. None when the files are the collected set.
.gs_seed_problems <- function(manifest, frozen) {
  if (file.exists(GLOW_SEED_MARKER))
    return(paste0("a collection of the pre-built files began and did not end (", GLOW_SEED_MARKER, " is there)"))
  stated <- unlist(frozen$artifacts)
  if (!length(stated)) return("the settings file states no checksums of the pre-built files (`artifacts`)")
  present <- .gs_seed_files(manifest)
  both <- intersect(names(stated), present)
  now <- unname(tools::md5sum(file.path(manifest$source, both)))
  other <- both[is.na(now) | now != unname(stated[both])]
  as.character(c(
    if (length(setdiff(names(stated), present))) paste("missing:", paste(setdiff(names(stated), present), collapse = ", ")),
    if (length(setdiff(present, names(stated)))) paste("not of the collected set:", paste(setdiff(present, names(stated)), collapse = ", ")),
    if (length(other)) paste("not the file that was collected:", paste(other, collapse = ", "))))
}

# The settings a seeded part stands for: its own and those of every part before it.
#
# parts  The table of parts.
# id     The part, as step.part.
#
# Returns the names: settings as section.name, and `content of <setting>` for an
# outside file and for an input path that a user may set in place of a part's
# output. A seed is valid for a run only when none of them differs. The content of
# an outside file is a checksum over every byte (.gs_setting_text()), so a raw input
# that was changed in place makes the supplied files invalid, and the same input
# under another path does not.
.gs_seed_closure <- function(parts, id, seen = character(0)) {
  spec <- parts[[id]]
  paths <- c(vapply(spec$files, `[[`, character(1), "setting"), unlist(lapply(spec$inputs, `[[`, "setting")))
  keys <- c(spec$settings, unlist(spec$region_settings), spec$forwarded, if (length(paths)) paste("content of", paths))
  for (inp in spec$inputs)
    if (!inp$part %in% c(seen, id)) keys <- c(keys, .gs_seed_closure(parts, inp$part, c(seen, id)))
  unique(keys)
}

# The content of data-example/seed-settings.yml: for each seeded part, the settings
# its files were built with (as text) and the content of its outside files (as
# checksums; "none" for an input path that was left to the layout).
#
# config  The configuration the files were built with. Default: the shipped example.
# cfg     The same configuration, already read (with absolute paths). The collector of
#         the pre-built files passes the object it validated the units with, so the
#         settings it writes cannot be those of a file that was edited meanwhile.
#
# Returns a list named by part, each list(settings, inputs). Called by
# data-example/write-seed-settings.R when the files are built, and by the test that
# keeps the shipped file in step with the shipped example.
.gs_seed_settings <- function(config = "example.yml", cfg = NULL) {
  if (is.null(cfg))
    cfg <- GLOWpipeline::glow_absolute_paths(GLOWpipeline::read_glow_config(
      config, base = "base-config/base_config.yml", registry = "base-config/settings-registry.yml", root = "."))
  parts <- GLOWpipeline::glow_read_parts("base-config/parts.yml", "base-config/settings-registry.yml")
  run <- .gs_run_new(cfg, GLOWpipeline::glow_layout(cfg), parts, NULL)
  manifest <- yaml::read_yaml("data-example/seed-manifest.yml")
  out <- list()
  for (id in unique(vapply(manifest$files, `[[`, character(1), "part"))) {
    keys <- sort(.gs_seed_closure(parts, id), method = "radix")     # the same order in every locale
    is_input <- startsWith(keys, "content of ")
    out[[id]] <- list(settings = as.list(.gs_setting_text(run, keys[!is_input])),
                      inputs = as.list(.gs_setting_text(run, keys[is_input])))
  }
  out
}

# One option of a stage script of release 0.3.6 that is a setting.
#
# The command line of a stage script names the configuration and the unit of work.
# It carries no setting: what a run does is in the configuration, and so it is in
# the completion record. A few values that release 0.3.6 took as options
# (--max-genes, --max-windows, --top-k, --n-cores, --allow-incomplete, --drilldown,
# --method) are settings of the configuration, and passing one of them stops the
# script with the name of the setting.
#
# st       The context that glow_stage_begin() returned.
# option   The name of the option, without the dashes.
# setting  c(section, name) of the setting that holds the value.
# default  The value when the setting gives none.
#
# Returns the setting's value, or `default` when it is null.
glow_stage_option <- function(st, option, setting = NULL, default = NULL) {
  given <- st$opts[[option]]
  is_given <- !is.null(given) && !identical(given, FALSE)
  if (is_given) .gs_refuse_option(option, if (is.null(setting)) NULL else paste(setting, collapse = "."))
  value <- if (is.null(setting)) NULL else st$cfg[[setting[1L]]][[setting[2L]]]
  if (is.null(value)) default else value
}

# The folders that the workflow keeps in a step folder for itself.
.gs_reserved_names <- c("records", "jobs", "logs", "provenance")

# Names that would make outputs collide.
#
# The registry type `file_name` keeps an output inside its folder. Two more things have
# to hold, which depend on what a name is for. A file-name pattern that a part writes
# once per chromosome has to contain {chr}: without it every chromosome writes the same
# file, and the last one stays. And the folder of the PI models, which is emptied before
# the training writes, must not be a folder that the workflow keeps for itself.
#
# cfg  The configuration.
.gs_check_names <- function(cfg) {
  for (key in c("gds_pattern", "favor_gds_pattern", "favor_csv_pattern")) {
    value <- cfg$data_prep[[key]]
    if (!is.null(value) && !grepl("{chr}", value, fixed = TRUE))
      stop("`data_prep.", key, "` is `", value, "`. It names one file per chromosome, so it has to contain {chr}: ",
           "without it every chromosome would write the same file.", call. = FALSE)
  }
  sub <- cfg$training$pi_models_subdir
  if (!is.null(sub) && sub %in% .gs_reserved_names)
    stop("`training.pi_models_subdir` is `", sub, "`, a folder that the workflow keeps for itself (",
         paste(.gs_reserved_names, collapse = ", "), "). Give the PI models a folder of their own.", call. = FALSE)
  invisible(NULL)
}

# Whether the run has not started in a folder: no completion record, and no unit that
# holds its lock. (A released lock of an earlier attempt is not a start.)
.gs_folder_untouched <- function(folder) {
  rec <- file.path(folder, "records")
  if (length(list.files(rec, pattern = "\\.yml$"))) return(FALSE)
  locks <- list.files(rec, pattern = "\\.lock\\.[0-9]+$", full.names = TRUE)
  !any(file.exists(file.path(locks, "owner")))
}

# Stop for an option of release 0.3.6 that was passed with a YAML configuration.
#
# option   The option, without the dashes.
# setting  The setting that holds the value. Default: the one this table names.
.gs_refuse_option <- function(option, setting = NULL) {
  known <- c("max-genes" = "snv_set.max_genes", "max-windows" = "snv_set.max_windows", "top-k" = "snv_set.top_k",
             "n-cores" = "shared.n_cores", "allow-incomplete" = "snv_set.allow_incomplete",
             "drilldown" = "snv_set.drilldown", "method" = "single_variant.calibration_method")
  if (is.null(setting) && option %in% names(known)) setting <- known[[option]]
  stop("--", option, " is not an option of a YAML configuration. ",
       if (!is.null(setting)) paste0("Set `", setting, "` in the configuration. ") else "",
       "The command line names the configuration and the unit of work, and carries no setting.", call. = FALSE)
}

#################### THE RECORDS OF A RUN ####################

# The state of one run: the configuration as the user wrote it (paths absolute,
# derived inputs still null), the layout, the table of parts, and what was looked
# up so far.
.gs_run_new <- function(cfg, lay, parts, region_type) {
  run <- new.env()
  run$cfg <- cfg; run$layout <- lay; run$parts <- parts; run$region_type <- region_type
  run$located <- list(); run$fingerprints <- list(); run$ancestors <- NULL; run$notes <- character(0)
  run$version <- if (file.exists(".release-source")) readLines(".release-source", n = 1L, warn = FALSE) else "development tree"
  run
}

# The value of a setting by its name in the table of parts: `dataset`,
# `section.name`, or `section.name.field` for a field of a map.
.gs_value <- function(cfg, key) {
  k <- strsplit(key, ".", fixed = TRUE)[[1]]
  value <- cfg
  for (name in k) {
    if (!is.list(value)) return(NULL)
    value <- value[[name]]
  }
  value
}

# A configuration with one setting replaced, by the same names.
.gs_set_value <- function(cfg, key, value) {
  k <- strsplit(key, ".", fixed = TRUE)[[1]]
  if (length(k) == 2L) cfg[[k[1]]][k[2]] <- list(value)
  else if (length(k) == 3L) cfg[[k[1]]][[k[2]]][k[3]] <- list(value)
  else stop("A derived setting is `section.name` or `section.name.field`: ", key)
  cfg
}

# The analyses up the `extends` chain that share the dataset and the output root,
# nearest first. Each is list(analysis, layout). The files of the chain are read
# once per run.
.gs_ancestors <- function(run) {
  if (!is.null(run$ancestors)) return(run$ancestors)
  files <- attr(run$cfg, "glow_chain")$file          # the base first, the user's file last
  out <- list()
  if (length(files) > 2L) for (f in rev(files[-c(1L, length(files))])) {
    anc <- GLOWpipeline::glow_absolute_paths(GLOWpipeline::read_glow_config(
      f, base = "base-config/base_config.yml", registry = "base-config/settings-registry.yml", root = "."))
    same <- identical(anc$dataset, run$cfg$dataset) && identical(anc$output_root, run$cfg$output_root)
    known <- c(run$cfg$analysis, vapply(out, `[[`, character(1), "analysis"))
    if (same && !anc$analysis %in% known)
      out[[length(out) + 1L]] <- list(analysis = anc$analysis,
                                      layout = GLOWpipeline::glow_layout(anc, region_types = run$cfg$snv_set$region_type))
  }
  run$ancestors <- out
  out
}

# The folder that holds the records of a part, in a layout.
.gs_part_folder <- function(lay, spec, region_type) {
  switch(spec$folder,
    "data.gds"             = lay$data$gds,
    "data.gds_favor"       = lay$data$gds_favor,
    "data.pcs"             = lay$data$pcs,
    "steps.data_prep"      = lay$steps$data_prep,
    "steps.training_b"     = lay$steps$training_b,
    "steps.training_pi"    = lay$steps$training_pi,
    "steps.single_variant" = lay$steps$single_variant,
    "steps.snv_set"        = file.path(lay$steps$snv_set_dir, region_type),
    stop("base-config/parts.yml: unknown folder `", spec$folder, "` of the part `", spec$id, "`."))
}

# The unit of work that a process is, from the selectors on its command line.
# A part that is split runs one unit per process, so its selector is required.
.gs_unit_name <- function(spec, opts, region_type) {
  need <- function(option, what) {
    if (is.null(opts[[option]]))
      stop("The part `", spec$part, "` of the step `", spec$step, "` runs once per ", what, ": pass --", option,
           ", as `glow run ", spec$step, "` does.", call. = FALSE)
    opts[[option]]
  }
  mode <- function() if (tolower(need("spa", "mode")) %in% c("on", "true", "1")) "spa" else "nospa"
  kind <- spec$units
  if (kind == "region") kind <- if (identical(region_type, "window")) "chunk" else "chromosome"
  switch(kind,
    none       = NULL,
    chromosome = paste0("chr", as.integer(need("chr", "chromosome"))),
    chunk      = sprintf("chunk%04d", as.integer(need("chunk", "chunk"))),
    group_mode = sprintf("%s.group%04d", mode(), as.integer(need("arrayid", "gene group"))),
    chunk_mode = sprintf("%s.chunk%04d", mode(), as.integer(need("chunk", "chunk"))))
}

# Whether the conditions of a `when` or a `run_when` of the table of parts hold.
#
# conds        The conditions (GLOWpipeline::glow_read_parts()): each names a setting
#              and what it has to be (`is`) or must not be (`is_not`): `null`, `set`
#              (not null and not empty) or the text of a value. Or it names the
#              region types it holds for.
# cfg          The configuration.
# region_type  The region type of the run, or NULL.
#
# Returns an empty string when all hold, else the first that does not, in words.
.gs_unmet <- function(conds, cfg, region_type) {
  for (cond in conds) {
    if (!is.null(cond$region_type)) {
      if (is.null(region_type) || !region_type %in% cond$region_type)
        return(sprintf("a part of the %s run only", paste(cond$region_type, collapse = " and ")))
      next
    }
    value <- .gs_value(cfg, cond$setting)
    is_set <- !is.null(value) && length(value) > 0L
    want <- if (is.null(cond$is)) cond$is_not else cond$is
    holds <- switch(want, null = !is_set, set = is_set,
                    is_set && is.atomic(value) && length(value) == 1L &&
                      identical(tolower(as.character(value)), tolower(want)))
    if (!is.null(cond$is_not)) holds <- !holds
    if (!holds) {
      shown <- if (!is_set) "null" else if (is.atomic(value) && length(value) == 1L) as.character(value) else "set"
      return(sprintf("`%s` is %s", cond$setting, shown))
    }
  }
  ""
}

# The units that a part has when it is finished. They follow from the configuration,
# never from the records that happen to exist.
#
# run     The state of the run.
# spec    The part.
# folder  The folder of the part: a window run reads its chunks from the table that
#         the preparation wrote there.
#
# Returns list(NULL) for a part that is not split, the names of the units for one
# that is, and character(0) when the units cannot be known yet (a window part
# before the preparation wrote the table of chunks).
.gs_expected_units <- function(run, spec, folder) {
  cfg <- run$cfg
  kind <- spec$units
  if (kind == "region") kind <- if (identical(run$region_type, "window")) "chunk" else "chromosome"
  chunks <- function() {
    f <- file.path(folder, "shared", "01-genome_chunks.rds")
    if (!file.exists(f)) return(character(0))
    sprintf("chunk%04d", as.integer(readRDS(f)$chunk_id))
  }
  # the native comparison runs every unit once per mode
  per_mode <- function(units) unlist(lapply(as.character(cfg$snv_set$staar_modes), function(m)
    if (length(units)) paste(m, units, sep = ".")))
  switch(kind,
    none       = list(NULL),
    chromosome = paste0("chr", cfg$shared$chroms),
    chunk      = chunks(),
    chunk_mode = as.character(per_mode(chunks())),
    group_mode = as.character(per_mode(sprintf("group%04d", .gs_gene_groups(cfg)))))
}

# The gene groups of the native coding comparison, as STAARpipeline numbers them:
# groups of `gene_num_in_array` genes per chromosome, counted over chromosomes 1 to
# 22. Returns the numbers of the groups of the run's chromosomes.
.gs_gene_groups <- function(cfg) {
  if (!requireNamespace("STAARpipeline", quietly = TRUE))
    stop("The native coding comparison needs the STAARpipeline package (its table of genes).", call. = FALSE)
  genes_info <- get("genes_info", envir = asNamespace("STAARpipeline"))
  groups <- ceiling(table(factor(genes_info[, 2], levels = 1:22)) / cfg$snv_set$gene_num_in_array)
  last <- cumsum(groups)
  as.integer(unlist(lapply(cfg$shared$chroms, function(k)
    if (groups[[k]] > 0) seq.int(last[[k]] - groups[[k]] + 1L, last[[k]]))))
}

# The arguments of a stage script that select a unit, from the unit's name.
.gs_unit_args <- function(unit) {
  if (is.null(unit)) return(character(0))
  if (grepl("^chr[0-9]+$", unit)) return(c("--chr", sub("^chr", "", unit)))
  if (grepl("^chunk[0-9]+$", unit)) return(c("--chunk", as.character(as.integer(sub("^chunk", "", unit)))))
  m <- regmatches(unit, regexec("^(spa|nospa)\\.(group|chunk)([0-9]+)$", unit))[[1]]
  if (length(m) != 4L) stop("Not the name of a unit: ", unit)
  c(if (m[3L] == "group") "--arrayid" else "--chunk", as.character(as.integer(m[4L])),
    "--spa", if (m[2L] == "spa") "on" else "off")
}

# The part of a run that reads the output of a part from the layout, or NULL.
#
# A part is left out of a run when its conditions (`run_when`) do not hold, which is
# how a part whose output the user supplies by a path is skipped. One supplied path
# need not replace every reader, though: with the annotated GDS files supplied, the
# principal components may still be computed from the plain GDS files, which the
# conversion from PLINK writes. So a part whose conditions do not hold still belongs
# to the run while a part of the run reads its output and no path of the user stands
# in for it.
#
# run   The state of the run.
# spec  The part whose conditions do not hold.
#
# Returns the name of a part of the same step that needs it, or NULL.
.gs_needed_by <- function(run, spec, seen = character(0)) {
  for (q in Filter(function(q) identical(q$step, spec$step) && !(q$id %in% c(seen, spec$id)), run$parts)) {
    reads <- Filter(function(r) identical(r$part, spec$id) && is.null(r$supplied), .gs_unit_reads(run, q, NULL))
    if (!length(reads)) next
    if (!nzchar(.gs_unmet(q$run_when, run$cfg, run$region_type)) || !is.null(.gs_needed_by(run, q, c(seen, spec$id))))
      return(q$part)
  }
  NULL
}

# The units of the producing part that an input reads, for a unit of the reading part.
#
# An input may name the one unit it reads (`unit`, with a setting in braces), or the
# setting that lists the chromosomes whose units it reads (`units_from`). Without
# either, a unit of one chromosome reads the unit of the same chromosome where the
# producing part has one, and every other unit reads all units.
#
# Returns the unit names, or NULL for all units of the producing part.
.gs_input_units <- function(run, inp, unit) {
  cfg <- run$cfg
  if (!is.null(inp$unit)) {
    wanted <- inp$unit
    for (k in regmatches(wanted, gregexpr("(?<=\\{)[^}]+(?=\\})", wanted, perl = TRUE))[[1]])
      wanted <- sub(paste0("{", k, "}"), as.character(.gs_value(cfg, k)), wanted, fixed = TRUE)
    return(wanted)
  }
  if (!is.null(inp$units_from)) return(paste0("chr", .gs_value(cfg, inp$units_from)))
  if (!is.null(unit) && run$parts[[inp$part]]$units %in% c("chromosome", "region") && grepl("^chr[0-9]+$", unit)) return(unit)
  NULL
}

# The outputs of other parts that a unit reads, as this run resolves them.
#
# One function says what a unit reads, for the completion records
# (.gs_unit_dependencies()) and for the SLURM mode (glow_slurm_plan() of
# lib/glow-slurm.R), so what a job waits for and what a starting unit checks cannot
# differ. An input of the table is dropped when its condition (`when`) does not hold
# in this run, or when its derived path does not apply.
#
# Returns a list, one element per input that is read: `part` (the producing part),
# `units` (the units of it that are read, NULL for all), `partial` (whether a unit
# that is not finished is tolerated), and `setting`, `pattern` and `supplied`: the
# path that the user set in place of the part's output, or NULL.
.gs_unit_reads <- function(run, spec, unit) {
  cfg <- run$cfg; rules <- .gs_derived_rules(); out <- list()
  for (inp in spec$inputs) {
    if (nzchar(.gs_unmet(inp$when, cfg, run$region_type))) next
    supplied <- if (!is.null(inp$setting)) .gs_value(cfg, inp$setting)
    if (is.null(supplied) && !is.null(inp$setting)) {
      rule <- rules[[inp$setting]]
      if (!is.null(rule$when) && !isTRUE(rule$when(cfg))) next
    }
    out[[length(out) + 1L]] <- list(part = inp$part, units = .gs_input_units(run, inp, unit), setting = inp$setting,
                                    pattern = inp$pattern, supplied = supplied,
                                    partial = !is.null(inp$partial) && isTRUE(.gs_value(cfg, inp$partial)))
  }
  out
}

# The real path of a file or folder that may not exist yet: the real path of the
# deepest folder that exists, with the rest of the path after it. Two spellings of
# one place (through a link, or with `..`) give the same text.
.gs_real_path <- function(path) {
  rest <- character(0)
  while (!file.exists(path) && dirname(path) != path) { rest <- c(basename(path), rest); path <- dirname(path) }
  do.call(file.path, as.list(c(normalizePath(path), rest)))
}

# The settings that the table names for a part in a run of a region type.
.gs_setting_keys <- function(spec, region_type)
  unique(c(spec$settings, spec$region_settings[[if (is.null(region_type)) "" else region_type]], spec$forwarded))

# Settings as the text that a record holds, named by setting. The value is the one
# the run uses (.gs_effective_value()). `content of <setting>` gives a checksum over
# every byte of the file or folder that the setting names, or "none": a seed record
# keeps the outside files this way, so it does not depend on where a clone lies or on
# the times of its files.
.gs_setting_text <- function(run, keys) {
  if (!length(keys)) return(stats::setNames(character(0), character(0)))
  vapply(keys, function(k) {
    if (startsWith(k, "content of ")) {
      v <- .gs_value(run$cfg, sub("^content of ", "", k))
      if (is.null(v)) return("none")
      key <- paste("content", v)
      if (is.null(run$fingerprints[[key]]))
        run$fingerprints[[key]] <- GLOWpipeline::glow_file_fingerprint(v, md5_limit = Inf, small = Inf)
      return(run$fingerprints[[key]])
    }
    GLOWpipeline::glow_setting_text(.gs_effective_value(run, k))
  }, character(1))
}

# The files that a unit reads from a folder the user named.
#
# A folder of per-chromosome files is identified by the files that are read, so a
# file added for another chromosome leaves a finished unit complete.
#
# run          The state of the run.
# path         The value of the setting: a file or a folder.
# pattern_key  The setting that holds the file-name pattern, or NULL.
# units        The units of the producing part that are read, as chr<N>, or NULL for
#              the chromosomes of the run.
#
# Returns `path` itself for a file, for a folder without a pattern and for a
# pattern without {chr}. Otherwise the files of the folder whose name is the
# pattern for one of the chromosomes, alone or followed by a dot and an extension
# (a PLINK file set is named without its extension). A chromosome without a file
# gives the name that is missing.
.gs_selected_files <- function(run, path, pattern_key, units = NULL) {
  if (is.null(pattern_key) || !dir.exists(path)) return(path)
  pattern <- .gs_effective_value(run, pattern_key)
  if (is.null(pattern) || !grepl("{chr}", pattern, fixed = TRUE)) return(path)
  chroms <- if (length(units) && all(grepl("^chr[0-9]+$", units))) sub("^chr", "", units) else run$cfg$shared$chroms
  wanted <- vapply(as.character(chroms), function(k) gsub("{chr}", k, pattern, fixed = TRUE), character(1))
  have <- list.files(path)
  of_name <- function(n) have[have == n | startsWith(have, paste0(n, "."))]
  hits <- lapply(wanted, of_name)
  file.path(path, unique(c(unlist(hits), wanted[lengths(hits) == 0L])))
}

# What a unit depends on: the four named lists of text that its completion record
# holds (GLOWpipeline::glow_unit_begin()).
#
# run   The state of the run.
# id    The part, as step.part.
# unit  The unit, or NULL.
# home  The analysis in which the part is looked at, as list(analysis, layout).
#       NULL: the run's own. A part of the same folder is read in that analysis.
# note  TRUE for the unit that runs: an input found up the chain is noted as reused.
#
# settings  every setting the table names for the part, as text, with the value the
#           run uses
# inputs    the generation mark of the record of each part that is read. Where the
#           parts have units of the same name, a unit depends on its own unit of
#           the other part. An input may name the chromosomes whose units it reads
#           (`units_from`). A part that no analysis holds finished gives "missing".
#           An input that may be partial (`partial`, a setting that is true) gives
#           "absent" for a unit that is not finished.
# files     the fingerprint of each outside file, and of each input path that the
#           user set in place of a part's output. From a folder with a file-name
#           pattern, the files that are read.
# versions  the packages and the workflow templates
#
# An input or outside file whose `when` does not hold is not read by this run and
# is left out. The list also carries `setting_text`, a function that gives the
# current text of any settings, for the comparison with a seed record.
.gs_unit_dependencies <- function(run, id, unit, home = NULL, note = FALSE) {
  spec <- run$parts[[id]]; cfg <- run$cfg
  fingerprint <- function(paths) {
    key <- paste(paths, collapse = "\n")
    if (is.null(run$fingerprints[[key]])) run$fingerprints[[key]] <- GLOWpipeline::glow_file_fingerprint(paths)
    run$fingerprints[[key]]
  }
  reads <- function(x) !nzchar(.gs_unmet(x$when, cfg, run$region_type))
  settings <- .gs_setting_text(run, .gs_setting_keys(spec, run$region_type))
  inputs <- character(0); files <- character(0)
  own_chromosome <- if (!is.null(unit) && grepl("^chr[0-9]+$", unit)) unit else NULL
  for (f in Filter(reads, spec$files)) {
    path <- .gs_value(cfg, f$setting)
    if (!is.null(path)) files[[f$setting]] <- fingerprint(.gs_selected_files(run, path, f$pattern, own_chromosome))
  }
  for (r in .gs_unit_reads(run, spec, unit)) {
    if (!is.null(r$supplied)) {                       # a file or folder of the user in place of the part
      files[[r$setting]] <- fingerprint(.gs_selected_files(run, r$supplied, r$pattern, r$units))
      next
    }
    found <- .gs_locate(run, r$part, units = r$units, from = spec, home = home, note = note, partial = r$partial)
    marks <- if (isTRUE(found$found)) found$generations else c("missing")
    names(marks) <- if (is.null(names(marks)) || !any(nzchar(names(marks)))) r$part else
      paste(r$part, names(marks), sep = ".")
    inputs <- c(inputs, marks)
  }
  list(settings = settings, inputs = inputs[!duplicated(names(inputs))], files = files,
       versions = GLOWpipeline::glow_run_versions(run$version),
       setting_text = function(keys) .gs_setting_text(run, keys))
}

# The analysis that holds a finished part.
#
# A run looks for the output of a part first in its own analysis and then in the
# analyses up its `extends` chain, nearest first. A candidate counts when every unit
# that is read has a record there, the record's outputs exist, and the record states
# the dependencies that the current configuration gives that part.
#
# run      The state of the run.
# id       The part that is looked for, as step.part.
# units    The units that are read. NULL: all units of the part.
# from     The part that reads it. A part of the same folder is looked for in one
#          analysis only (`home`), because the scripts of one folder read each other's
#          files in place. A part of the dataset level has one folder for all analyses.
# home     The analysis of the part that reads it. NULL: the run's own.
# note     Whether a part found up the chain is noted as reused.
# partial  TRUE: a unit that is not finished does not disqualify the candidate. Its
#          mark is "absent".
#
# Returns list(found, analysis, layout, folder, generations, why). `generations` is
# named by unit. `why` says what is missing or different in the nearest candidate.
.gs_locate <- function(run, id, units = NULL, from = NULL, home = NULL, note = FALSE, partial = FALSE) {
  spec <- run$parts[[id]]
  own <- list(analysis = run$cfg$analysis, layout = run$layout)
  same_folder <- !is.null(from) && identical(from$folder, spec$folder) && identical(from$level, spec$level)
  candidates <- if (spec$level == "dataset") list(own)
                else if (same_folder) list(if (is.null(home)) own else home)
                else c(list(own), .gs_ancestors(run))
  key <- paste(id, paste(units, collapse = ","), paste(vapply(candidates, `[[`, character(1), "analysis"), collapse = ">"),
               if (partial) "partial" else "")
  # A result that is already known is taken from the cache. The note of a reuse is made in
  # both cases: the first look at a part is often made without a note, while the dependencies
  # of another part are followed, and the record of the folder still has to state the reuse.
  noted <- function(result) {
    if (note && isTRUE(result$found) && !identical(result$analysis, run$cfg$analysis))
      run$notes <- union(run$notes, sprintf("%s/%s: reused from %s", spec$step, spec$part, result$analysis))
    result
  }
  if (!is.null(run$located[[key]])) return(noted(run$located[[key]]))
  why <- NULL; result <- NULL
  for (cand in candidates) {
    folder <- .gs_part_folder(cand$layout, spec, run$region_type)
    expected <- if (is.null(units)) .gs_expected_units(run, spec, folder) else as.list(units)
    if (!length(expected)) {
      if (is.null(why)) why <- sprintf("in the analysis `%s` its units are not known yet", cand$analysis)
      next
    }
    marks <- character(0); ok <- TRUE
    for (u in expected) {
      state <- GLOWpipeline::glow_unit_state(folder, spec$part, u, .gs_unit_dependencies(run, id, u, home = cand))
      if (!state$complete && partial) { marks[[if (is.null(u)) "" else u]] <- "absent"; next }
      if (!state$complete) {
        if (is.null(why) || !is.null(state$record))
          why <- sprintf("in the analysis `%s`%s %s%s", cand$analysis, if (is.null(u)) "" else paste0(", unit ", u, ","),
                         state$reason, if (length(state$differences)) paste0(":\n", paste0("      - ", state$differences, collapse = "\n")) else "")
        ok <- FALSE; break
      }
      marks[[if (is.null(u)) "" else u]] <- state$record$generation
    }
    if (ok) {
      result <- list(found = TRUE, analysis = cand$analysis, layout = cand$layout, folder = folder,
                     generations = marks, why = "")
      break
    }
  }
  if (is.null(result)) result <- list(found = FALSE, analysis = NULL, layout = NULL, folder = NULL,
                                      generations = character(0), why = if (is.null(why)) "it has no completion record" else why)
  run$located[[key]] <- result
  noted(result)
}

# Whether the folder of a part is finished in an analysis up the chain.
#
# The scripts of one folder read each other's files in place, so a folder is reused
# whole or not at all. It is reusable from an analysis when every part of the folder
# that belongs to this run (the part that asks, and each part whose `run_when`
# holds) is finished there in all its units, with the dependencies of the current
# configuration. What has to be finished follows from the configuration, never from
# the records that exist: a folder with a finished preparation and no scan is not a
# finished folder.
#
# Returns list(analysis, folder) of the nearest such analysis, or NULL.
.gs_folder_reusable <- function(run, spec) {
  of_run <- Filter(function(p) identical(p$folder, spec$folder) && identical(p$level, spec$level) &&
                     (identical(p$id, spec$id) || !nzchar(.gs_unmet(p$run_when, run$cfg, run$region_type))), run$parts)
  for (anc in .gs_ancestors(run)) {
    folder <- .gs_part_folder(anc$layout, spec, run$region_type)
    if (!length(list.files(file.path(folder, "records"), pattern = "\\.yml$"))) next
    ok <- TRUE
    for (p in of_run) {
      expected <- .gs_expected_units(run, p, folder)
      if (!length(expected)) { ok <- FALSE; break }
      for (u in expected)
        if (!GLOWpipeline::glow_unit_state(folder, p$part, u, .gs_unit_dependencies(run, p$id, u, home = anc))$complete) { ok <- FALSE; break }
      if (!ok) break
    }
    if (ok) return(list(analysis = anc$analysis, folder = folder))
  }
  NULL
}

# A folder that one unit owns must hold the files of its record, and no other.
#
# The SNV-set step reads every model of the PI models folder. That folder is the
# output of one unit, whose record lists its files. A file in the folder that the
# record does not list (a model copied in by hand, or left by another version of the
# run) would join the ensemble without any record of it. So when a part reads, by a
# derived path, a folder that is the output of a part with one unit, the run stops if
# the folder holds a file that the record does not list.
.gs_check_owned_folders <- function(run, spec) {
  cfg <- run$cfg
  for (inp in spec$inputs) {
    if (is.null(inp$setting) || !is.null(.gs_value(cfg, inp$setting))) next       # no path, or the user's own folder
    if (nzchar(.gs_unmet(inp$when, cfg, run$region_type))) next
    other <- run$parts[[inp$part]]
    if (!identical(other$units, "none")) next
    found <- .gs_locate(run, inp$part, from = spec)
    if (!isTRUE(found$found)) next
    dir <- .gs_derived_rules()[[inp$setting]]$path(found$layout, cfg, run)
    if (!dir.exists(dir)) next
    rec <- GLOWpipeline::glow_unit_record(found$folder, other$part, NULL)
    listed <- file.path(found$folder, rec$outputs[!grepl("^/", rec$outputs)])
    extra <- setdiff(list.files(dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE), listed)
    if (length(extra))
      stop("The folder ", dir, " holds ", length(extra), " file(s) that the completion record of `", other$step, "/", other$part,
           "` does not list: ", paste(utils::head(basename(extra), 5L), collapse = ", "), ". `", spec$step, "/", spec$part,
           "` would read them as results of that part. Remove them, or run `glow run ", other$step, " --part ", other$part,
           " --force` to produce the folder again.", call. = FALSE)
  }
  invisible(NULL)
}

# The two checks on the inputs of a unit, before any work.
#
# The missing input: a part whose output this unit reads, and that no analysis
# holds finished for this configuration, stops the run. The message names the step
# and the part to run first and says what is missing or different.
#
# The example-data check: when `dataset` is not `example`, an input that lies
# inside data-example/, or inside a dataset folder that carries the marker
# EXAMPLE-DATA, stops the run and names the setting. What is checked is what the
# unit reads: the file, the files that are read from a folder with a file-name
# pattern, or every file of a folder. Paths are compared in their real form, so a
# link into the example data is caught, also a link inside a folder. A model
# trained on the example would look normal in a real study.
.gs_check_inputs <- function(run, spec, unit, deps) {
  .gs_check_owned_folders(run, spec)
  missing <- names(deps$inputs)[deps$inputs == "missing"]
  if (length(missing)) {
    lines <- vapply(missing, function(m) {
      id <- names(run$parts)[vapply(names(run$parts), function(p) identical(m, p) || startsWith(m, paste0(p, ".")), logical(1))]
      id <- id[which.max(nchar(id))]
      p <- run$parts[[id]]
      found <- Filter(function(x) identical(x$found, FALSE), run$located[startsWith(names(run$located), paste0(id, " "))])
      found <- found[order(-nchar(vapply(found, `[[`, character(1), "why")))]     # the most telling reason first
      # A part of the unit's own folder that an analysis up the chain holds finished: this analysis
      # reuses that folder whole, so an ordinary run of the part reports the reuse and writes
      # nothing here. A unit that runs again in this analysis needs the part in this analysis.
      hint <- ""
      if (identical(p$folder, spec$folder) && identical(p$level, spec$level) && p$level != "dataset") {
        from <- .gs_folder_reusable(run, p)
        if (!is.null(from))
          hint <- sprintf(paste0("\n      The analysis `%s` holds this folder finished, and this analysis reuses it whole. To run a part ",
                                 "of it again here, give this analysis its own results of the parts before it: ",
                                 "`glow run %s --part %s --force`, or `glow run %s --force` for the whole step."),
                          from$analysis, p$step, p$part, p$step)
      }
      sprintf("  - %s, part %s (`glow run %s`): %s%s", p$step, p$part, p$step, if (length(found)) found[[1]]$why else "not finished", hint)
    }, character(1))
    stop("`", spec$step, "/", spec$part, "` reads the results of ", length(missing), " part(s) that are not finished for this ",
         "configuration:\n", paste(unique(lines), collapse = "\n"), "\nRun them first.", call. = FALSE)
  }
  if (identical(run$cfg$dataset, "example")) return(invisible(NULL))
  example_root <- normalizePath("data-example", mustWork = FALSE)
  known <- new.env()
  dir_in_example <- function(dir) {
    if (!is.null(known[[dir]])) return(known[[dir]])
    hit <- startsWith(dir, paste0(example_root, "/")) || identical(dir, example_root)
    up <- dir
    while (!hit && nchar(up) > 1L) {
      if (file.exists(file.path(up, "EXAMPLE-DATA"))) hit <- TRUE
      if (identical(dirname(up), up)) break
      up <- dirname(up)
    }
    known[[dir]] <- hit
    hit
  }
  in_example <- function(path) {
    real <- normalizePath(path, mustWork = FALSE)
    dirs <- if (dir.exists(real))
      unique(c(real, dirname(normalizePath(list.files(real, recursive = TRUE, full.names = TRUE), mustWork = FALSE))))
      else dirname(real)
    any(vapply(dirs, dir_in_example, logical(1)))
  }
  cfg <- .gs_derive_inputs(run, spec)
  reads <- function(x) !nzchar(.gs_unmet(x$when, run$cfg, run$region_type))
  own_chromosome <- if (!is.null(unit) && grepl("^chr[0-9]+$", unit)) unit else NULL
  consumed <- list()
  for (x in Filter(reads, spec$files)) {
    value <- .gs_value(cfg, x$setting)
    if (is.character(value) && length(value) == 1L) consumed[[x$setting]] <- .gs_selected_files(run, value, x$pattern, own_chromosome)
  }
  for (x in Filter(function(i) !is.null(i$setting) && reads(i), spec$inputs)) {
    value <- .gs_value(cfg, x$setting)
    if (is.character(value) && length(value) == 1L)
      consumed[[x$setting]] <- .gs_selected_files(run, value, x$pattern, .gs_input_units(run, x, unit))
  }
  bad <- names(consumed)[vapply(consumed, function(paths) any(vapply(paths, in_example, logical(1))), logical(1))]
  if (length(bad))
    stop("The dataset is `", run$cfg$dataset, "`, and ", length(bad), " input(s) of `", spec$step, "/", spec$part,
         "` belong to the shipped example:\n",
         paste0("  - ", bad, ": ", vapply(bad, function(k) .gs_value(cfg, k), character(1)), collapse = "\n"),
         "\nSet these to the files of your own dataset. The example's files are valid for `dataset: example` only.",
         call. = FALSE)
  invisible(NULL)
}

# config-resolved.yml of a step folder: every setting with the value the parts of
# the folder use, every derived path in absolute form, the reuse notes and the
# provenance. A derived input that no part of the folder reads stays null.
# The file is written to a temporary folder and renamed, so units that run at the
# same time never leave half a file.
.gs_write_resolved <- function(run, cfg, dir, what) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = "resolved-", tmpdir = dir); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE, expand = FALSE), add = TRUE)
  GLOWpipeline::write_resolved_config(cfg, tmp, provenance = c(
    list(glowanalyses = run$version, written_by = what,
         command = paste(c("Rscript", sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)),
                           commandArgs(trailingOnly = TRUE)), collapse = " ")),
    if (length(run$notes)) list(reuse = run$notes)))
  if (!file.rename(file.path(tmp, "config-resolved.yml"), file.path(dir, "config-resolved.yml")))
    warning("Could not write ", file.path(dir, "config-resolved.yml"), call. = FALSE)
  invisible(file.path(dir, "config-resolved.yml"))
}

# End a script whose unit has nothing to do: say why, and quit with status 0.
# With options(glow.stage.no_quit = TRUE), which the tests set, the context is
# returned with `skipped` instead.
.gs_finish_early <- function(st, reason, message) {
  cat(message, "\n", sep = "")
  st$skipped <- reason
  if (!isTRUE(getOption("glow.stage.no_quit"))) quit(save = "no", status = 0L)
  st
}
