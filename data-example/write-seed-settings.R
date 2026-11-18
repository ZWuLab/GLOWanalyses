#!/usr/bin/env Rscript
# ==============================================================================
# File: data-example/write-seed-settings.R
# Purpose: Write data-example/seed-settings.yml, the settings that the example's
#   pre-built files (data-example/seed-manifest.yml) were built with. The file is
#   written when those files are built and is not edited by hand: it fixes what the
#   files stand for, so that a later edit of example.yml does not change it. On a
#   run of the example, lib/glow-stage.R writes these settings into the seed
#   records, and a run whose settings differ is told to produce its own files.
#
# For each seeded part the file holds
#   settings  every setting of the part and of the parts before it, as the text
#             that a completion record holds (GLOWpipeline::glow_setting_text())
#   inputs    the content of the outside files of those parts, as a checksum over
#             every byte, and of the input paths that a user may set in place of a
#             part's output ("none" when the path is derived). A raw input of the
#             example that is changed in place therefore makes the pre-built files
#             invalid, and a clone at another path does not.
#
# Usage (from the GLOWanalyses directory):
#   Rscript data-example/write-seed-settings.R [--collect | --check]
#
#   --collect  first copy the files of data-example/seed-manifest.yml from the run
#              of example.yml back into data-example/cohort/, then write the file.
#              Every unit of the manifest has to be finished under example.yml as it
#              is now, by a run and not by supplied files: otherwise nothing is
#              copied and the script stops. A folder is replaced by the files that
#              its completion record lists.
#   --check    write nothing: compare the file with what example.yml and the pre-built
#              files give now, and exit with status 1 when they differ or when a
#              collection was interrupted.
#
# The file also holds `artifacts`, the checksum of every pre-built file. A run of the
# example uses the files only when they are the set that the file states
# (.gs_seed_problems() of lib/glow-stage.R). While --collect puts files in place, the
# marker data-example/seed-collection.incomplete is there. It is removed after the file
# is written, so a collection that was interrupted leaves it, and neither --check nor a
# run of the example accepts the files until a collection has run to its end.
#
# Run it only through data-example/build-intermediates.sh, which rebuilds the files
# from example.yml and calls this script with --collect.
source("lib/glow-stage.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(args, c("--check", "--collect"))))
  stop("Unknown argument(s): ", paste(setdiff(args, c("--check", "--collect")), collapse = " "),
       ". Usage: Rscript data-example/write-seed-settings.R [--collect | --check]", call. = FALSE)
check <- "--check" %in% args
collect <- "--collect" %in% args
# --check promises to write nothing, and --collect replaces the tracked pre-built files.
if (check && collect) stop("Use --collect or --check, not both: --check writes nothing.", call. = FALSE)
out_file <- "data-example/seed-settings.yml"
built <- NULL              # the settings of the pre-built files
artifacts <- NULL          # the checksums of the pre-built files, named by their path under the cohort
publish <- function() invisible(NULL)
manifest <- yaml::read_yaml("data-example/seed-manifest.yml")

# ---- --collect: the rebuilt files, from the run of example.yml into the tracked cohort ----
# The files are copied into a staging folder first. They enter the tracked cohort only when,
# after the copy, every unit is still finished by the record it was validated with, under
# dependencies that are looked up anew, and neither a file of the configuration nor an input
# of the files has changed: a unit or a part before it that is run again meanwhile, a raw
# input that is replaced meanwhile, or an example.yml that is edited meanwhile, publishes
# nothing. The settings that the file states are computed before the copy, so they are those
# of the inputs that the validated units were made from.
if (collect) {
  watched <- c("example.yml", "base-config/base_config.yml", "base-config/settings-registry.yml", "base-config/parts.yml",
               "data-example/seed-manifest.yml")
  sums_before <- tools::md5sum(watched)
  cfg <- GLOWpipeline::glow_absolute_paths(GLOWpipeline::read_glow_config(
    "example.yml", base = "base-config/base_config.yml", registry = "base-config/settings-registry.yml", root = "."))
  lay <- GLOWpipeline::glow_layout(cfg)
  parts <- GLOWpipeline::glow_read_parts("base-config/parts.yml", "base-config/settings-registry.yml")
  run <- .gs_run_new(cfg, lay, parts, NULL)
  # the files that were read are the files that were checksummed, and they stay so to the end
  unchanged <- function(when) {
    chain <- attr(cfg, "glow_chain")
    now <- tools::md5sum(c(watched, chain$file))
    was <- c(sums_before, stats::setNames(chain$md5, chain$file))
    changed <- unique(names(now)[is.na(now) | unname(now) != unname(was)])
    if (length(changed))
      stop("A file of the configuration changed ", when, ": ", paste(changed, collapse = ", "),
           ". Nothing was collected. Run data-example/build-intermediates.sh again.", call. = FALSE)
  }
  unchanged("while it was read")
  if (!identical(cfg$analysis, manifest$analysis))
    stop("example.yml names the analysis `", cfg$analysis, "`, and data-example/seed-manifest.yml the analysis `",
         manifest$analysis, "`. The files of one cannot be collected as the other.", call. = FALSE)
  # every unit of the manifest: finished under this configuration, by a run
  # `run` keeps what it looked up (the checksums of the outside files, the records of the parts
  # that a unit reads), so the state after the copy is asked of a new one
  state_of <- function(item, run) {
    spec <- parts[[item$part]]
    GLOWpipeline::glow_unit_state(.gs_part_folder(lay, spec, NULL), spec$part, item$unit, .gs_unit_dependencies(run, item$part, item$unit))
  }
  listed <- list()
  for (item in manifest$files) {
    folder <- .gs_part_folder(lay, parts[[item$part]], NULL)
    state <- state_of(item, run)
    what <- paste0(item$part, if (!is.null(item$unit)) paste0(" ", item$unit))
    if (!state$complete)
      stop("The unit ", what, " is not finished under example.yml as it is now (", state$reason,
           "). Nothing was copied. Run data-example/build-intermediates.sh to its end.", call. = FALSE)
    if (isTRUE(state$record$seeded))
      stop("The unit ", what, " holds supplied files, not the result of a run. Nothing was copied. ",
           "data-example/build-intermediates.sh runs the two steps with --force.", call. = FALSE)
    listed[[length(listed) + 1L]] <- list(item = item, what = what, generation = state$record$generation,
                                          files = file.path(folder, state$record$outputs[!grepl("^/", state$record$outputs)]))
  }
  # the settings and the checksums of the inputs, as they are while the units are valid
  built <- .gs_seed_settings("example.yml", cfg = cfg)
  # the staging folder is in the temporary folder of this R process, which R removes at its end,
  # also when the script stops: nothing is left in the tree
  stage <- tempfile(pattern = "collect-"); dir.create(stage)
  for (l in listed) {
    from <- file.path(lay$dataset_dir, l$item$to); staged <- file.path(stage, l$item$from)
    if (dir.exists(from)) {
      members <- l$files[startsWith(l$files, paste0(from, "/"))]        # the files of the record, no other
      rel <- substring(members, nchar(from) + 2L)
      for (d in unique(dirname(rel))) dir.create(file.path(staged, d), recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(members, file.path(staged, rel), overwrite = TRUE)
    } else {
      if (!from %in% l$files) stop("The completion record of ", l$item$part, " does not list ", from, ".", call. = FALSE)
      dir.create(dirname(staged), recursive = TRUE, showWarnings = FALSE)
      ok <- file.copy(from, staged, overwrite = TRUE)
    }
    if (!all(ok)) stop("Could not copy ", from, " to ", staged, ".", call. = FALSE)
  }
  # After the copy: the same records, valid under dependencies that are looked up now, so the
  # bytes that were copied are those of units that are finished for the inputs as they are.
  fresh <- .gs_run_new(cfg, lay, parts, NULL)
  for (l in listed) {
    state <- state_of(l$item, fresh)
    if (!state$complete || isTRUE(state$record$seeded) || !identical(state$record$generation, l$generation))
      stop("The unit ", l$what, " was run again, lost its completion record, or is no longer finished for its inputs (",
           if (state$complete) "another record" else state$reason, "), while its files were copied. ",
           "Nothing was collected. Run data-example/build-intermediates.sh again, and nothing else at the same time.", call. = FALSE)
  }
  if (!identical(.gs_seed_settings("example.yml", cfg = cfg), built))
    stop("An input of the pre-built files changed while they were copied. Nothing was collected. ",
         "Run data-example/build-intermediates.sh again, and nothing else at the same time.", call. = FALSE)
  unchanged("while the files were copied")
  staged_files <- list.files(stage, recursive = TRUE, all.files = TRUE, no.. = TRUE)
  artifacts <- stats::setNames(unname(tools::md5sum(file.path(stage, staged_files))), staged_files)
  # The staged files enter the tracked cohort after the settings file is ready to be written
  # (below). The marker is there from the first replacement until the settings file is written.
  publish <- function() {
    writeLines(c("A collection of the example's pre-built files began and has not ended.",
                 "While this file is there, the files under data-example/cohort/ and data-example/seed-settings.yml",
                 "may belong to two collections. Run data-example/build-intermediates.sh to its end."), GLOW_SEED_MARKER)
    for (l in listed) {
      to <- file.path(manifest$source, l$item$from); staged <- file.path(stage, l$item$from)
      unlink(to, recursive = TRUE, expand = FALSE)
      dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
      ok <- if (dir.exists(staged)) file.copy(staged, dirname(to), recursive = TRUE) else file.copy(staged, to, overwrite = TRUE)
      if (!all(ok)) stop("Could not put ", staged, " in place as ", to, ".", call. = FALSE)
      cat(sprintf("collected %s -> %s\n", file.path(lay$dataset_dir, l$item$to), to))
    }
  }
}

# ---- without --collect: the settings of example.yml as it is, and the files in place ----
if (!collect) {
  if (!check && file.exists(GLOW_SEED_MARKER))
    stop(GLOW_SEED_MARKER, " is there: a collection of the pre-built files began and did not end. ",
         "Run data-example/build-intermediates.sh to its end, which collects them again.", call. = FALSE)
  built <- .gs_seed_settings("example.yml")
  in_place <- .gs_seed_files(manifest)
  artifacts <- stats::setNames(unname(tools::md5sum(file.path(manifest$source, in_place))), in_place)
}
artifacts <- artifacts[order(names(artifacts), method = "radix")]      # the same order in every locale
quote_text <- function(x) paste0("\"", gsub("\"", "\\\\\"", gsub("\\\\", "\\\\\\\\", x)), "\"")
lines <- c(
  "# ==============================================================================",
  "# seed-settings.yml - the settings the example's pre-built files were built with",
  "#",
  "# Written by data-example/write-seed-settings.R when the files of",
  "# data-example/seed-manifest.yml are built. Do not edit: this file fixes what the",
  "# pre-built files stand for. Each value is the text of a setting as a completion",
  "# record holds it. `inputs` are the outside files, by a checksum of their content.",
  "# `artifacts` are the pre-built files themselves, by path under the cohort folder and",
  "# checksum: a run of the example uses them only when they are this set.",
  "# ==============================================================================")
for (id in names(built)) {
  lines <- c(lines, paste0(id, ":"), "  settings:")
  for (k in names(built[[id]]$settings))
    lines <- c(lines, sprintf("    %s: %s", k, quote_text(built[[id]]$settings[[k]])))
  lines <- c(lines, "  inputs:")
  for (k in names(built[[id]]$inputs))
    lines <- c(lines, sprintf("    %s: %s", quote_text(k), quote_text(built[[id]]$inputs[[k]])))
}

lines <- c(lines, "artifacts:", sprintf("  %s: %s", quote_text(names(artifacts)), quote_text(unname(artifacts))))

if (check) {
  interrupted <- file.exists(GLOW_SEED_MARKER)
  same <- file.exists(out_file) && identical(readLines(out_file, warn = FALSE), lines)
  ok <- same && !interrupted
  cat(sprintf("Verdict for the seed settings: %s\n", if (ok) "PASS" else "FAIL"))
  if (interrupted) cat(GLOW_SEED_MARKER, "is there: a collection of the pre-built files began and did not end.\n")
  if (!same) cat(out_file, "does not state the settings that example.yml gives now, or the pre-built files that are in place.",
                 "Rebuild the pre-built files (data-example/build-intermediates.sh), which writes it.\n")
  quit(save = "no", status = if (ok) 0L else 1L)
}
if (collect) { unchanged("before the files were put in place"); publish() }
tmp <- paste0(out_file, ".tmp-", Sys.getpid())
writeLines(lines, tmp)
if (!file.rename(tmp, out_file)) { unlink(tmp, expand = FALSE); stop("Could not write ", out_file, ".", call. = FALSE) }
# the files and the settings file are one collection again
if (collect) unlink(GLOW_SEED_MARKER, expand = FALSE)
cat(sprintf("Wrote %s: %d seeded part(s), %d pre-built file(s).\n", out_file, length(built), length(artifacts)))
