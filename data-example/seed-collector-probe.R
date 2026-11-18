#!/usr/bin/env Rscript
# ==============================================================================
# File: data-example/seed-collector-probe.R
# Purpose: Run data-example/write-seed-settings.R --collect with one disturbance, as a
#   second process could cause it while the collector works. It is the helper of
#   data-example/check-seed-collector.sh, which runs it in a copy of the GLOWanalyses
#   tree whose builder has finished, and never in the tree itself: the disturbances
#   change files of the example and of its run.
#
# Usage (from the copy of the GLOWanalyses directory):
#   Rscript data-example/seed-collector-probe.R <kind>
#
#   edit       example.yml is edited at the first copy into the staging folder
#   rerun      the B model is produced again at that moment: its record is gone and
#              its file holds other bytes
#   raw        a raw input of the example (the phenotype table) is changed at that
#              moment
#   upstream   a part that the pre-built files were made from, and that is not itself
#              collected (the principal components), is produced again at that moment:
#              its record has another generation
#   interrupt  the collector stops after it has put its first file in place
#
# The disturbed files are left as they are. The calling script puts them back.
.probe_kind <- commandArgs(trailingOnly = TRUE)[1]   # (the collector has a variable `what` of its own)
if (is.na(.probe_kind) || !.probe_kind %in% c("edit", "rerun", "raw", "upstream", "interrupt"))
  stop("Usage: Rscript data-example/seed-collector-probe.R edit|rerun|raw|upstream|interrupt", call. = FALSE)
.probe_done <- FALSE
.probe_run <- "runs/example"                         # the run of example.yml in the copy

# One disturbance, at the first copy into the staging folder.
.probe_disturb <- function() {
  if (.probe_done) return(invisible()); .probe_done <<- TRUE
  if (.probe_kind == "edit") {
    x <- readLines("example.yml"); writeLines(sub("pi_n_models: 5", "pi_n_models: 6", x, fixed = TRUE), "example.yml")
    cat("PROBE: example.yml edited (pi_n_models 5 -> 6) during the copy\n")
  } else if (.probe_kind == "rerun") {
    rec <- file.path(.probe_run, "baseline/training/b/records/b_estimate.yml")
    file.rename(rec, paste0(rec, ".kept"))
    file.copy(file.path(.probe_run, "baseline/training/b/b_model.rds"), "b_model.kept", overwrite = TRUE)
    writeLines("INCOMPLETE PRODUCER OUTPUT", file.path(.probe_run, "baseline/training/b/b_model.rds"))
    cat("PROBE: the record of b_estimate removed and its model overwritten during the copy\n")
  } else if (.probe_kind == "raw") {
    cat("\n", file = "data-example/cohort/pheno.csv", append = TRUE)
    cat("PROBE: data-example/cohort/pheno.csv changed during the copy\n")
  } else if (.probe_kind == "upstream") {
    rec <- list.files(.probe_run, pattern = "^compute_pcs.*[.]yml$", recursive = TRUE, full.names = TRUE)
    rec <- unique(normalizePath(rec[grepl("/records/", rec)]))    # the analysis folder links to the data of the dataset
    if (length(rec) != 1L) stop("PROBE: expected one record of compute_pcs under ", .probe_run, ", found ", length(rec))
    x <- readLines(rec)
    hit <- grep("^\"generation\":", x)
    if (length(hit) != 1L) stop("PROBE: the record ", rec, " has no generation line")
    x[hit] <- "\"generation\": \"produced-again-by-the-probe\""
    writeLines(x, rec)
    cat("PROBE: the record of compute_pcs got another generation during the copy (", rec, ")\n", sep = "")
  }
}
if (.probe_kind == "interrupt") {
  # after the first file has entered the tracked cohort
  trace(file.copy, print = FALSE, where = baseenv(), exit = quote(
    if (any(startsWith(as.character(to), "data-example/cohort"))) stop("PROBE: interrupted after the first file was put in place", call. = FALSE)))
} else {
  trace(file.copy, tracer = quote(.probe_disturb()), print = FALSE, where = baseenv())
}
commandArgs <- function(trailingOnly = FALSE) "--collect"
source("data-example/write-seed-settings.R")
