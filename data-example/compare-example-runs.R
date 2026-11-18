#!/usr/bin/env Rscript
# ==============================================================================
# File: data-example/compare-example-runs.R
# Purpose: Compare two runs of the shipped example, file by file and cell by cell.
#   It answers one question: did a change of the workflow (a new release, the YAML
#   form of the configuration, the `glow` command, the Terra tasks) leave the
#   results of the example as they were?
#
#   Each run is a GLOWanalyses directory in which the example was run, in one of
#   two forms: "R" (the config.R files of release 0.3.6, outputs beside each
#   config.R) or "yaml" (one YAML configuration, outputs on the derived layout
#   runs/example/). The script maps the folders of the two forms onto each other,
#   lists the files of both runs, and compares
#     - every table (.csv, .tsv): the same columns, the same rows, and every cell
#       identical as text, so a missing cell has to be missing in both;
#     - every R object that is a result (.rds): identical, or equal without any
#       numeric difference when it holds a function;
#     - every GDS file: the same nodes with identical data;
#     - the two run records of an SNV-set run (01-config_used.rds,
#       01-input_checksums.rds) field by field.
#   Values that describe the run and not the data (a creation time, a run name, a
#   package version, a path) cannot be equal between two runs. They are listed in
#   RUN_KEYS, RUN_COLUMNS and OBJECT_RUN_FIELDS below, left out of the comparison,
#   and, for the keys, checked on run B instead: a run name has to be the name of
#   run B, a path has to exist, and the GDS folder and the phenotype bundle have to
#   be those of run B.
#
#   A PASS needs more than no difference. Each group that is asked for has to hold
#   the results that the example produces (REQUIRED below), in both runs, so two
#   runs that both stopped early do not pass. A checksum record has to state the
#   checksum of the files it names. And the version of GLOWr that a file carries is
#   left out only when it is the version that the run itself recorded.
#
# Usage (from any directory):
#   Rscript data-example/compare-example-runs.R \
#     --a <GLOWanalyses dir> --a-form R|yaml --b <GLOWanalyses dir> --b-form R|yaml \
#     [--groups single_variant,gene,window,coding,training,data_prep] [--out <file.md>] \
#     [--tolerance <relative tolerance>] [--a-root <output root>] [--b-root <output root>]
#
#   --tolerance  Two runs on one machine give the same text in every cell, and the
#                default asks for that (0). Two runs on different machines can differ
#                in the last digits of a number, because the linear algebra library
#                takes other code paths on other processors. A tolerance such as 1e-8
#                compares numbers by their relative difference and everything else as
#                text. The report says which rule was used and the largest difference.
#   --a-root, --b-root  The output root of a YAML run, when it is not runs/ in the
#                GLOWanalyses directory.
#
# Output: a Markdown report (to --out, or the screen) and, as the last line on the
# screen, "Verdict for the comparison: PASS" or "... FAIL". Exit status 1 on FAIL.

suppressMessages(library(data.table))

#################### WHAT IS COMPARED, AND WHAT IS NOT ####################

# The folders of one run, by group. `form` is "R" or "yaml".
run_folders <- function(dir, form, root = NULL) {
  if (form == "R") {
    ex <- file.path(dir, "runs", "example"); co <- file.path(dir, "data-example", "cohort")
    list(single_variant = file.path(ex, "single-variant", "outputs"),
         ld_scores      = file.path(co, "example_ld"),
         gene           = file.path(ex, "snv-set", "outputs"),
         window         = file.path(ex, "snv-set-window", "outputs"),
         coding         = file.path(ex, "snv-set-coding", "outputs"),
         training       = file.path(ex, "training", "outputs"),
         gds            = file.path(co, "example_gds"),
         gds_favor_gds  = file.path(co, "example_gds_favor", "flexible", "gds"),
         gds_favor_csv  = file.path(co, "example_gds_favor", "flexible", "csv"),
         pcs            = file.path(co, "example_pcs"),
         pheno          = file.path(co, "example_pheno"))
  } else {
    ds <- file.path(if (is.null(root)) file.path(dir, "runs") else root, "example"); an <- file.path(ds, "baseline")
    list(single_variant = file.path(an, "single_variant"),
         ld_scores      = file.path(an, "single_variant", "ld_scores"),
         gene           = file.path(an, "snv_set", "gene"),
         window         = file.path(an, "snv_set", "window"),
         coding         = file.path(an, "snv_set", "coding"),
         training       = c(file.path(an, "training", "b"), file.path(an, "training", "pi")),
         gds            = file.path(ds, "data", "gds"),
         gds_favor_gds  = file.path(ds, "data", "gds_favor", "flexible", "gds"),
         gds_favor_csv  = file.path(ds, "data", "gds_favor", "flexible", "csv"),
         pcs            = file.path(ds, "data", "pcs"),
         pheno          = file.path(an, "data_prep"))
  }
}

# The groups behind each name that --groups accepts.
GROUPS <- list(single_variant = c("single_variant", "ld_scores"), gene = "gene", window = "window",
               coding = "coding", training = "training",
               data_prep = c("gds", "gds_favor_gds", "gds_favor_csv", "pcs", "pheno"))

# What a finished run of the example holds, per group: each pattern has to match a
# file that both runs have and that was compared. The example has chromosomes 21
# and 22. Without this, two runs that both lack a result would not differ.
REQUIRED <- list(
  single_variant = c("results/single_variant_chr21\\.csv$", "results/single_variant_chr22\\.csv$",
                     "results/single_variant_all\\.csv$", "results/single_variant_calibration\\.tsv$"),
  ld_scores      = "(^|/)ld_scores\\.rds$",
  gene           = c("results/glow_chr21\\.csv$", "results/glow_chr22\\.csv$", "shared/01-config_used\\.rds$",
                     "shared/01-input_checksums\\.rds$", "aggregated/.+\\.csv$"),
  window         = c("results/scan_chunk_[0-9]+\\.csv$", "shared/01-genome_chunks\\.rds$", "shared/01-config_used\\.rds$",
                     "aggregated/scan_results_all\\.csv$"),
  coding         = c("results-glow/glow_coding_chr21\\.csv$", "results-glow/glow_coding_chr22\\.csv$",
                     "shared/01-config_used\\.rds$", "aggregated/glow_coding_all\\.csv$"),
  # the metadata of the PI training states how many models it trained (see stated_results())
  training       = c("(^|/)b_model\\.rds$", "pi_models/model_1\\.rds$", "pi_models/training_metadata\\.rds$"),
  gds            = c("chr21.*\\.gds$", "chr22.*\\.gds$"),
  gds_favor_gds  = c("chr21.*\\.gds$", "chr22.*\\.gds$"),
  gds_favor_csv  = c("chr21.*\\.csv$", "chr22.*\\.csv$"),
  pcs            = "(^|/)pcs\\.rds$",
  pheno          = "(^|/)pheno_covar\\.rds$")

# Files that are records of the run, not results: logs, provenance, the completion
# records, the resolved configuration and the job scripts of a YAML run, plots. They are counted and
# not compared. A plot is a drawing of a table that is compared.
SKIP_PATTERN <- "(^|/)(logs|slurm-logs|provenance|records|jobs)/|(^|/)config-resolved\\.yml$|README\\.md$|\\.pdf$|\\.txt$|\\.log$"

# Results that exist only when a part ran one unit per chromosome: the LD scores of
# each chromosome, which the combining part puts together into ld_scores.rds. They
# are compared when both runs have them. In one run only, they are counted.
PER_UNIT_PATTERN <- "(^|/)ld_chr[0-9]+\\.rds$"

# One file name for both forms: the single-variant tables were named marginal_* in
# release 0.3.6, and the phenotype bundle carried the name of its run folder.
canonical_name <- function(rel) {
  rel <- sub("(^|/)marginal_", "\\1single_variant_", rel)
  rel <- sub("(^|/)data-prep_pheno_covar\\.rds$", "\\1pheno_covar.rds", rel)
  rel
}

# Keys of the key/value tables, and fields of the run records, that describe the run.
# They are left out of the comparison and checked on run B (check_run_value()).
RUN_KEYS <- c(
  "run_name", "created_at", "report_created_at", "prepared_at",           # the name and the time of the run
  "GFisher", "GLOWr", "GLOWpipeline", "GLOWanalyses", "git_sha", "session_info",  # the software
  "config_path", "config_source",                                          # the configuration file
  "output_dir", "gds_dir", "pheno_path", "pi_model_dir", "b_model_path",   # folders that the layout moves
  "scan_table", "ld_scores", "ld_scores_path", "calibrated_table", "report",
  "staar_null_path", "staar_null_spa_path", "staar_null_nospa_path", "staar_catalog_path")

# Columns of an ordinary table that describe the run.
RUN_COLUMNS <- c("run")

# Named parts of an R object that describe the run: the label of the run in the
# phenotype bundle, and folders and files in the records of the PCs and PI models.
OBJECT_RUN_FIELDS <- c("covar_set", "output_dir", "case_csv", "control_source", "gds_files", "timestamp")

# Attributes of a GDS node that describe the run (a folder on the machine).
GDS_RUN_ATTRIBUTES <- c("favor_db_path")

# The field and the GDS attribute in which GLOWr writes its own version into a file
# it produces (the B model, the annotated GDS files). Two runs with two versions of
# GLOWr differ there. The value is left out only when it is the GLOWr version that
# the run itself recorded (RUN_VERSION, from 01-run_info.tsv), so a file from another
# run does not pass as this run's.
VERSION_FIELD <- "glowr_version"

# Fields of 01-config_used.rds that the two forms write differently by design. They
# are compared by a rule, not by identity (compare_config_used()):
#   filter_spec  the YAML form writes all seven fields of the filter, the config.R
#                form the four it sets, so the built filters are compared
#   top_k        one example.yml sets top_k for all region types, and only the gene
#                summary reads it, so a difference counts for a gene run only

#################### READING AND COMPARING ####################

list_files <- function(folders) {
  out <- list()
  for (folder in folders) {
    if (!dir.exists(folder)) next
    rel <- list.files(folder, recursive = TRUE, all.files = TRUE, no.. = TRUE)
    # the LD scores have a group of their own
    rel <- rel[!grepl("^ld_scores/", rel)]
    # a link to a data folder of the dataset is not a file of this folder
    top <- unique(sub("/.*$", "", rel[grepl("/", rel)]))
    links <- top[nzchar(Sys.readlink(file.path(folder, top)))]
    if (length(links)) rel <- rel[!sub("/.*$", "", rel) %in% links]
    for (r in rel) out[[canonical_name(r)]] <- file.path(folder, r)
  }
  out
}

read_text_table <- function(path) {
  # every cell as the text that is in the file; an empty cell is missing
  as.data.frame(fread(path, colClasses = "character", na.strings = "", header = TRUE,
                      keepLeadingZeros = TRUE, showProgress = FALSE), stringsAsFactors = FALSE)
}

is_key_value <- function(d) identical(names(d), c("key", "value"))

# Compare two tables. Returns list(ok, rows, cols, note, run_values).
compare_tables <- function(pa, pb) {
  a <- read_text_table(pa); b <- read_text_table(pb)
  run_values <- NULL
  if (is_key_value(a) && is_key_value(b)) {
    keys_a <- a$key; keys_b <- b$key
    run_values <- data.frame(key = intersect(RUN_KEYS, union(keys_a, keys_b)), stringsAsFactors = FALSE)
    run_values$a <- a$value[match(run_values$key, keys_a)]
    run_values$b <- b$value[match(run_values$key, keys_b)]
    a <- a[!a$key %in% RUN_KEYS, , drop = FALSE]; b <- b[!b$key %in% RUN_KEYS, , drop = FALSE]
    # the checksum of the phenotype bundle, when the two bundles differ in their run label only
    if (BUNDLE$label_only && "pheno_bundle_md5" %in% a$key && "pheno_bundle_md5" %in% b$key &&
        identical(a$value[a$key == "pheno_bundle_md5"], BUNDLE$md5_a) &&
        identical(b$value[b$key == "pheno_bundle_md5"], BUNDLE$md5_b)) {
      a <- a[a$key != "pheno_bundle_md5", , drop = FALSE]; b <- b[b$key != "pheno_bundle_md5", , drop = FALSE]
      run_values <- rbind(run_values, data.frame(key = "pheno_bundle_md5", a = BUNDLE$md5_a, b = BUNDLE$md5_b,
                                                 stringsAsFactors = FALSE))
    }
  } else {
    drop <- intersect(RUN_COLUMNS, intersect(names(a), names(b)))
    if (length(drop)) {
      run_values <- data.frame(key = drop, a = vapply(drop, function(k) paste(unique(a[[k]]), collapse = " "), character(1)),
                               b = vapply(drop, function(k) paste(unique(b[[k]]), collapse = " "), character(1)),
                               stringsAsFactors = FALSE)
      a <- a[, setdiff(names(a), drop), drop = FALSE]; b <- b[, setdiff(names(b), drop), drop = FALSE]
    }
  }
  if (!identical(names(a), names(b)))
    return(list(ok = FALSE, rows = nrow(a), cols = ncol(a), run_values = run_values,
                note = sprintf("columns differ: only in A {%s}; only in B {%s}",
                               paste(setdiff(names(a), names(b)), collapse = ", "),
                               paste(setdiff(names(b), names(a)), collapse = ", "))))
  if (nrow(a) != nrow(b))
    return(list(ok = FALSE, rows = nrow(a), cols = ncol(a), run_values = run_values,
                note = sprintf("%d rows in A, %d in B", nrow(a), nrow(b))))
  largest <- 0
  differing <- vapply(names(a), function(col) {
    x <- a[[col]]; y <- b[[col]]
    differs <- xor(is.na(x), is.na(y)) | (!is.na(x) & !is.na(y) & x != y)
    if (TOLERANCE > 0 && any(differs)) {
      # a pair of numbers within the tolerance is not a difference
      xn <- suppressWarnings(as.numeric(x)); yn <- suppressWarnings(as.numeric(y))
      both <- differs & !is.na(xn) & !is.na(yn)
      rel <- abs(xn - yn) / pmax(abs(xn), abs(yn), .Machine$double.xmin)
      rel[!is.finite(rel)] <- ifelse(xn[!is.finite(rel)] == yn[!is.finite(rel)], 0, Inf)
      largest <<- max(largest, rel[both], 0, na.rm = TRUE)
      differs[both & rel <= TOLERANCE] <- FALSE
    }
    sum(differs)
  }, numeric(1))
  bad <- differing[differing > 0]
  note <- if (length(bad) == 0L) "" else
    paste(sprintf("%s: %d cell(s)", names(bad), bad), collapse = "; ")
  if (largest > 0) note <- paste0(note, if (nzchar(note)) "; " else "", sprintf("largest relative difference of a number %.2e", largest))
  list(ok = length(bad) == 0L, rows = nrow(a), cols = ncol(a), note = note, run_values = run_values,
       within = length(bad) == 0L && largest > 0)
}

# An R object without the named parts that describe the run, at any depth and in
# its attributes. Returns list(object, removed).
strip_run_fields <- function(x, side) {
  removed <- character(0)
  own_version <- function(v) !is.na(RUN_VERSION[[side]]) && identical(as.character(v), RUN_VERSION[[side]])
  walk <- function(v) {
    at <- attributes(v)
    for (k in setdiff(names(at), c("names", "class", "row.names", "dim", "dimnames", "levels"))) {
      if (is.list(at[[k]])) attr(v, k) <- walk(at[[k]])
    }
    if (is.list(v) && !is.null(names(v))) {
      # the version rule looks at the elements with that name only: another element can be
      # anything, a function of a fitted model included
      hit <- names(v) %in% OBJECT_RUN_FIELDS
      named_version <- which(names(v) == VERSION_FIELD)
      hit[named_version] <- hit[named_version] | vapply(v[named_version], own_version, logical(1))
      if (any(hit)) { removed <<- union(removed, names(v)[hit]); v <- v[!hit] }
      for (i in seq_along(v)) if (is.list(v[[i]])) v[i] <- list(walk(v[[i]]))
    }
    v
  }
  list(object = walk(x), removed = removed)
}

# Whether two R objects are equal within the tolerance, by the rule of the tables: every
# pair of numbers is compared on its own, by its relative difference, and everything that
# is not a number has to be equal without a numeric difference. (all.equal() with a
# tolerance compares the mean relative difference of a whole vector, so one wrong number
# among many that agree would pass.) Returns list(ok, largest): `largest` is the largest
# relative difference of a pair of numbers that was looked at.
objects_within <- function(a, b) {
  largest <- 0
  other_attributes <- function(v) { at <- attributes(v); at[setdiff(names(at), "names")] }
  walk <- function(x, y) {
    if (is.numeric(x) && is.numeric(y)) {
      if (length(x) != length(y) || !isTRUE(all.equal(attributes(x), attributes(y), tolerance = 0))) return(FALSE)
      if (!identical(is.na(x), is.na(y)) || !identical(is.nan(x), is.nan(y))) return(FALSE)
      here <- !is.na(x)
      rel <- abs(x[here] - y[here]) / pmax(abs(x[here]), abs(y[here]), .Machine$double.xmin)
      odd <- !is.finite(rel)                                   # two infinite values
      rel[odd] <- ifelse(x[here][odd] == y[here][odd], 0, Inf)
      largest <<- max(largest, rel, 0)
      return(all(rel <= TOLERANCE))
    }
    if (is.list(x) && is.list(y)) {
      if (length(x) != length(y) || !identical(names(x), names(y))) return(FALSE)
      if (!isTRUE(all.equal(other_attributes(x), other_attributes(y), tolerance = 0))) return(FALSE)
      for (i in seq_along(x)) if (!walk(x[[i]], y[[i]])) return(FALSE)
      return(TRUE)
    }
    isTRUE(all.equal(x, y, tolerance = 0))
  }
  list(ok = walk(a, b), largest = largest)
}

# Compare two R objects that are results. Verdicts: identical, or equal without any
# numeric difference (an object that holds a function is never identical to a second
# copy of itself), or equal within the tolerance (objects_within()), each possibly after
# leaving out the parts that describe the run.
compare_objects <- function(pa, pb) {
  a <- readRDS(pa); b <- readRDS(pb)
  if (identical(a, b)) return(list(ok = TRUE, verdict = "identical", note = ""))
  sa <- strip_run_fields(a, "a"); sb <- strip_run_fields(b, "b")
  left_out <- union(sa$removed, sb$removed)
  note <- if (length(left_out)) paste("left out as run values:", paste(left_out, collapse = ", ")) else ""
  if (identical(sa$object, sb$object)) return(list(ok = TRUE, verdict = "identical", note = note))
  eq <- all.equal(sa$object, sb$object, tolerance = 0)
  if (isTRUE(eq)) return(list(ok = TRUE, verdict = "equal, no numeric difference", note = note))
  if (TOLERANCE > 0) {
    within <- objects_within(sa$object, sb$object)
    if (within$ok)
      return(list(ok = TRUE, verdict = "equal within the tolerance",
                  note = paste(c(if (nzchar(note)) note, sprintf("largest relative difference of a number %.2e", within$largest)), collapse = "; ")))
  }
  list(ok = FALSE, verdict = "FAIL: differs", note = paste(c(note, utils::head(eq, 3L)), collapse = "; "))
}

# Compare two GDS files by content: the same nodes, each with identical data, and
# the same attributes apart from those that name a folder of the run.
compare_gds <- function(pa, pb) {
  if (identical(unname(tools::md5sum(pa)), unname(tools::md5sum(pb))))
    return(list(ok = TRUE, verdict = "identical (md5)", note = ""))
  if (!requireNamespace("gdsfmt", quietly = TRUE))
    return(list(ok = FALSE, verdict = "FAIL: differs (md5)", note = "gdsfmt is not installed, so the content was not compared"))
  a <- gdsfmt::openfn.gds(pa); on.exit(gdsfmt::closefn.gds(a), add = TRUE)
  b <- gdsfmt::openfn.gds(pb); on.exit(gdsfmt::closefn.gds(b), add = TRUE)
  nodes <- function(g, dirs) gdsfmt::ls.gdsn(g, recursive = TRUE, include.dirs = dirs, include.hidden = TRUE)
  na <- nodes(a, FALSE); nb <- nodes(b, FALSE)
  if (!identical(na, nb)) return(list(ok = FALSE, verdict = "FAIL: differs", note = "the node lists differ"))
  bad <- na[!vapply(na, function(n) identical(gdsfmt::read.gdsn(gdsfmt::index.gdsn(a, n)),
                                              gdsfmt::read.gdsn(gdsfmt::index.gdsn(b, n))), logical(1))]
  attr_bad <- character(0); left_out <- character(0)
  for (n in c("", nodes(a, TRUE))) {
    xa <- gdsfmt::get.attr.gdsn(if (n == "") a$root else gdsfmt::index.gdsn(a, n))
    xb <- gdsfmt::get.attr.gdsn(if (n == "") b$root else gdsfmt::index.gdsn(b, n))
    for (k in union(names(xa), names(xb))) if (!identical(xa[[k]], xb[[k]])) {
      own_versions <- k == VERSION_FIELD && !anyNA(unlist(RUN_VERSION)) &&
        identical(as.character(xa[[k]]), RUN_VERSION$a) && identical(as.character(xb[[k]]), RUN_VERSION$b)
      if (k %in% GDS_RUN_ATTRIBUTES || own_versions) left_out <- union(left_out, k) else attr_bad <- c(attr_bad, paste0(n, "@", k))
    }
  }
  ok <- length(bad) == 0L && length(attr_bad) == 0L
  note <- paste(c(sprintf("%d nodes", length(na)),
                  if (length(bad)) paste("data differ in:", paste(utils::head(bad, 5L), collapse = ", ")),
                  if (length(attr_bad)) paste("attributes differ:", paste(utils::head(attr_bad, 5L), collapse = ", ")),
                  if (length(left_out)) paste("left out as run values:", paste(left_out, collapse = ", "))), collapse = "; ")
  list(ok = ok, verdict = if (ok) "identical content" else "FAIL: differs", note = note)
}

# Compare the settings record of an SNV-set run field by field.
compare_config_used <- function(pa, pb) {
  a <- readRDS(pa); b <- readRDS(pb)
  fields <- union(names(a), names(b))
  run_values <- data.frame(key = intersect(RUN_KEYS, fields), stringsAsFactors = FALSE)
  whole <- function(x) if (is.null(x)) NA_character_ else paste(format(x), collapse = " ")   # shortened for the report later
  run_values$a <- vapply(run_values$key, function(k) whole(a[[k]]), character(1))
  run_values$b <- vapply(run_values$key, function(k) whole(b[[k]]), character(1))
  problems <- character(0); by_rule <- character(0)
  for (f in setdiff(fields, RUN_KEYS)) {
    if (identical(a[[f]], b[[f]])) next
    if (f == "filter_spec" && requireNamespace("GLOWr", quietly = TRUE) &&
        identical(do.call(GLOWr::variant_filter, a[[f]]), do.call(GLOWr::variant_filter, b[[f]]))) {
      by_rule <- c(by_rule, f); next
    }
    if (f == "top_k" && !identical(a$region_type, "gene")) { by_rule <- c(by_rule, f); next }
    problems <- c(problems, f)
  }
  note <- paste(c(if (length(problems)) paste("fields differ:", paste(problems, collapse = ", ")),
                  if (length(by_rule)) paste("equal by a declared rule:", paste(by_rule, collapse = ", "))),
                collapse = "; ")
  list(ok = length(problems) == 0L, note = note, run_values = run_values)
}

# Compare the checksum record of an SNV-set run: the same inputs by content.
# The record holds one checksum per input file (the phenotype bundle, the B model,
# the PI models, the GDS files). First each record is checked against its own run:
# every file it names has to exist and to have the checksum it states. Then the
# inputs are matched by file name. Returns the names whose checksums differ, with the
# two checksums, which the main part resolves: an input that the two runs built
# separately can differ in a run value only (a label, a time, a folder), and it then
# has a row of its own in this report that says so.
compare_checksums <- function(pa, pb, dir_a, dir_b) {
  # What a record has to list: the inputs that the settings record beside it names. A record
  # that lists other files, or not all of them, says nothing about the inputs of the run,
  # whatever its checksums are. Returns "" or what is wrong.
  inventory <- function(p, dir) {
    snap <- file.path(dirname(p), "01-config_used.rds")
    if (!file.exists(snap)) return("the settings record 01-config_used.rds is not beside it")
    s <- readRDS(snap); rec <- readRDS(p)
    real <- function(x) normalizePath(ifelse(startsWith(x, "/"), x, file.path(dir, x)), mustWork = FALSE)
    want <- list(
      pheno_path     = s$pheno_path,
      b_model_path   = s$b_model_path,
      pi_model_files = list.files(real(s$pi_model_dir), full.names = TRUE, pattern = "\\.rds$"),
      gds_files      = file.path(s$gds_dir, vapply(as.character(s$chroms), function(k) gsub("{chr}", k, s$gds_pattern, fixed = TRUE), character(1))))
    if (!is.list(rec) || !setequal(names(rec), names(want)))
      return(sprintf("it holds the entries %s, and the record of a run holds %s",
                     paste(names(rec), collapse = ", "), paste(names(want), collapse = ", ")))
    for (kind in names(want)) {
      have <- names(rec[[kind]]); have <- if (is.null(have)) character(0) else have[nzchar(have)]
      need <- if (is.null(want[[kind]])) character(0) else want[[kind]]
      extra <- have[!real(have) %in% real(need)]; lacking <- need[!real(need) %in% real(have)]
      if (anyDuplicated(real(have))) return(sprintf("its entry %s names a file twice", kind))
      if (length(extra)) return(sprintf("its entry %s names %s, which the run does not read by its settings record", kind, extra[1L]))
      if (length(lacking)) return(sprintf("its entry %s does not name %s, which the run reads by its settings record", kind, lacking[1L]))
    }
    ""
  }
  read <- function(p, dir) {
    v <- unlist(lapply(unname(readRDS(p)), function(e) stats::setNames(unname(e), names(e))))
    v <- v[!(is.na(v) & !nzchar(names(v)))]                      # an input that the run does not have (no B model file)
    if (!length(v)) return(list(sums = character(0), wrong = character(0)))
    paths <- ifelse(startsWith(names(v), "/"), names(v), file.path(dir, names(v)))
    actual <- unname(tools::md5sum(paths))                       # NA for a file that does not exist
    wrong <- names(v)[is.na(actual) | actual != v]
    sums <- stats::setNames(unname(v), basename(names(v)))
    list(sums = sums[order(names(sums))], wrong = wrong)
  }
  for (side in list(list("A", pa, dir_a), list("B", pb, dir_b))) {
    wrong <- inventory(side[[2]], side[[3]])
    if (nzchar(wrong))
      return(list(ok = FALSE, differing = character(0),
                  note = sprintf("the record of run %s does not list the inputs of its run: %s", side[[1]], wrong)))
  }
  ra <- read(pa, dir_a); rb <- read(pb, dir_b)
  for (side in list(list("A", ra), list("B", rb))) if (length(side[[2]]$wrong))
    return(list(ok = FALSE, differing = character(0),
                note = sprintf("the record of run %s does not state the checksum of the file it names: %s", side[[1]],
                               paste(utils::head(side[[2]]$wrong, 3L), collapse = ", "))))
  va <- ra$sums; vb <- rb$sums
  if (length(va) == 0L || !identical(names(va), names(vb)))
    return(list(ok = FALSE, differing = character(0),
                note = sprintf("the input lists differ: %d input(s) in A, %d in B", length(va), length(vb))))
  differing <- names(va)[va != vb]
  list(ok = length(differing) == 0L, differing = differing, sums_a = va[differing], sums_b = vb[differing],
       note = sprintf("%d input(s), each with the checksum its record states, %d with a different checksum in the two runs",
                      length(va), length(differing)))
}

# A value that describes the run, checked on run B.
#
# key, value  The key and the value in run B.
# group       The group of the file.
# b           What is known about run B: form, dir (its GLOWanalyses directory),
#             folders (run_folders()) and bundle (the phenotype bundle it reads).
#
# Returns what the value is, in words. The words in capitals fail the comparison.
check_run_value <- function(key, value, group, b) {
  same_path <- function(x, y) identical(normalizePath(x, mustWork = FALSE), normalizePath(y, mustWork = FALSE))
  in_run <- function(p) if (startsWith(p, "/")) p else file.path(b$dir, p)
  # the name of run B: its dataset and analysis, or the folder of its config.R
  run_name <- if (b$form == "yaml") "example/baseline" else {
    folder <- b$folders[[group]][1L]
    if (identical(basename(folder), "outputs")) basename(dirname(folder)) else NA_character_
  }
  if (key %in% c(RUN_COLUMNS, "run_name")) {
    if (is.na(value) || !nzchar(value)) return("MISSING")
    if (is.na(run_name)) return("a run name")
    return(if (identical(value, run_name)) "the name of run B" else "UNEXPECTED")
  }
  if (is.na(value) || !nzchar(value) || identical(value, "NA")) return(if (key %in% c("b_model_path", "git_sha", "config_source",
      "staar_null_path", "staar_null_spa_path", "staar_null_nospa_path", "staar_catalog_path", "report",
      "report_created_at", "ld_scores", "ld_scores_path")) "empty, allowed" else "MISSING")
  if (key == "pheno_bundle_md5") return("the checksum of run B's bundle")
  if (key %in% c("created_at", "report_created_at", "prepared_at"))
    return(if (!is.na(as.POSIXct(substr(value, 1L, 19L), tz = "UTC", optional = TRUE))) "a time" else "NOT A TIME")
  # the two inputs that every step reads have to be those of run B
  if (key == "gds_dir") return(if (same_path(in_run(value), b$folders$gds_favor_gds)) "the GDS folder of run B" else "NOT THE GDS FOLDER OF RUN B")
  if (key == "pheno_path") return(if (same_path(in_run(value), b$bundle)) "the bundle of run B" else "NOT THE BUNDLE OF RUN B")
  if (key %in% c("config_path", "output_dir", "pi_model_dir", "b_model_path", "scan_table", "ld_scores", "ld_scores_path",
                 "calibrated_table", "report", "staar_null_path", "staar_null_spa_path", "staar_null_nospa_path",
                 "staar_catalog_path"))
    # some tables keep the name of a file or folder, not its path: a name is the same in every run
    return(if (!grepl("/", value, fixed = TRUE)) "a name, not a path" else if (file.exists(in_run(value))) "an existing path" else "PATH NOT FOUND")
  "present"
}

#################### MAIN ####################

args <- commandArgs(trailingOnly = TRUE)
opt <- list(groups = paste(names(GROUPS)[1:4], collapse = ","), out = NA_character_, tolerance = "0")
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!startsWith(args[i], "--") || i == length(args)) stop("Unknown or dangling argument: ", args[i])
  opt[[key]] <- args[i + 1L]; i <- i + 2L
}
for (need in c("a", "a-form", "b", "b-form")) if (is.null(opt[[need]])) stop("Missing --", need)
stopifnot(opt[["a-form"]] %in% c("R", "yaml"), opt[["b-form"]] %in% c("R", "yaml"),
          dir.exists(opt$a), dir.exists(opt$b))
TOLERANCE <- suppressWarnings(as.numeric(opt$tolerance))
if (is.na(TOLERANCE) || TOLERANCE < 0 || TOLERANCE >= 1) stop("--tolerance is a relative tolerance, a number from 0 to below 1.")
asked <- strsplit(opt$groups, ",")[[1]]
if (!all(asked %in% names(GROUPS))) stop("--groups takes: ", paste(names(GROUPS), collapse = ", "))

fa <- run_folders(opt$a, opt[["a-form"]], opt[["a-root"]]); fb <- run_folders(opt$b, opt[["b-form"]], opt[["b-root"]])

# The phenotype bundle that the analysis steps of each run read. Its checksum is in
# the records of those steps. Two bundles that are identical apart from the label of
# their run (`covar_set`) have different checksums, which is then not a difference
# of the results.
bundle_path <- function(dir, form, root = NULL) if (form == "R") file.path(dir, "data-example", "cohort", "pheno_covar.rds") else
  file.path(if (is.null(root)) file.path(dir, "runs") else root, "example", "baseline", "data_prep", "pheno_covar.rds")
# The GLOWr version that each run recorded for itself, in the run table of an SNV-set
# preparation (01-run_info.tsv, key GLOWr). NA when the run has none.
RUN_VERSION <- lapply(list(a = fa, b = fb), function(folders) {
  for (f in file.path(unlist(folders[c("gene", "window", "coding")]), "shared", "01-run_info.tsv")) if (file.exists(f)) {
    info <- read_text_table(f)
    if ("GLOWr" %in% info$key) return(info$value[info$key == "GLOWr"][1L])
  }
  NA_character_
})
BUNDLE <- list(label_only = FALSE, md5_a = NA_character_, md5_b = NA_character_)
local({
  pa <- bundle_path(opt$a, opt[["a-form"]], opt[["a-root"]]); pb <- bundle_path(opt$b, opt[["b-form"]], opt[["b-root"]])
  if (file.exists(pa) && file.exists(pb)) {
    BUNDLE$md5_a <<- unname(tools::md5sum(pa)); BUNDLE$md5_b <<- unname(tools::md5sum(pb))
    if (!identical(BUNDLE$md5_a, BUNDLE$md5_b)) {
      a <- readRDS(pa); b <- readRDS(pb)
      BUNDLE$label_only <<- !identical(a$covar_set, b$covar_set) &&
        identical(a[names(a) != "covar_set"], b[names(b) != "covar_set"])
    }
  }
})
rows <- list(); run_rows <- list(); pending_checksums <- list(); paths <- list()
for (group in unlist(GROUPS[asked], use.names = FALSE)) {
  la <- list_files(fa[[group]]); lb <- list_files(fb[[group]])
  if (group == "single_variant") {   # in the R form the LD scores sit in the data folder
    la <- la[!grepl("^ld_", basename(names(la)))]; lb <- lb[!grepl("^ld_", basename(names(lb)))]
  }
  for (name in sort(union(names(la), names(lb)))) {
    pa <- la[[name]]; pb <- lb[[name]]
    row <- list(group = group, file = name, kind = "", size = "", verdict = "", note = "")
    if (grepl(SKIP_PATTERN, name)) {
      row$kind <- "record or plot"
      row$verdict <- if (!is.null(pa) && !is.null(pb)) "not compared, in both" else
        sprintf("not compared, only in %s", if (is.null(pb)) "A" else "B")
    } else if ((is.null(pa) || is.null(pb)) && grepl(PER_UNIT_PATTERN, name)) {
      row$kind <- "per-chromosome result"
      row$verdict <- sprintf("not compared, only in %s", if (is.null(pb)) "A" else "B")
      row$note <- "the run that has it ran this part one chromosome at a time; the combined result is compared"
    } else if (is.null(pa) || is.null(pb)) {
      row$kind <- "result"; row$verdict <- sprintf("FAIL: only in %s", if (is.null(pb)) "A" else "B")
    } else if (grepl("\\.(csv|tsv)$", name)) {
      res <- compare_tables(pa, pb)
      row$kind <- "table"; row$size <- sprintf("%d x %d", res$rows, res$cols)
      row$verdict <- if (!res$ok) "FAIL: differs" else if (isTRUE(res$within)) "equal within the tolerance" else "identical"
      row$note <- res$note
      if (!is.null(res$run_values) && nrow(res$run_values)) run_rows[[length(run_rows) + 1L]] <- cbind(group = group, file = name, res$run_values)
    } else if (grepl("01-config_used\\.rds$", name)) {
      res <- compare_config_used(pa, pb)
      row$kind <- "settings record"; row$verdict <- if (res$ok) "equal" else "FAIL: differs"; row$note <- res$note
      if (nrow(res$run_values)) run_rows[[length(run_rows) + 1L]] <- cbind(group = group, file = name, res$run_values)
    } else if (grepl("01-input_checksums\\.rds$", name)) {
      res <- compare_checksums(pa, pb, opt$a, opt$b)
      row$kind <- "checksum record"; row$verdict <- if (res$ok) "equal" else "FAIL: differs"; row$note <- res$note
      if (length(res$differing)) pending_checksums[[length(pending_checksums) + 1L]] <-
        list(row = length(rows) + 1L, differing = res$differing, sums_a = res$sums_a, sums_b = res$sums_b)
    } else if (grepl("\\.rds$", name)) {
      res <- compare_objects(pa, pb)
      row$kind <- "R object"; row$verdict <- res$verdict; row$note <- res$note
    } else if (grepl("\\.gds$", name)) {
      res <- compare_gds(pa, pb)
      row$kind <- "GDS file"; row$verdict <- res$verdict; row$note <- res$note
    } else {
      row$kind <- "other"; row$verdict <- "not compared, in both"
    }
    rows[[length(rows) + 1L]] <- as.data.frame(row, stringsAsFactors = FALSE)
    paths[[length(rows)]] <- c(a = if (is.null(pa)) NA_character_ else pa, b = if (is.null(pb)) NA_character_ else pb)
  }
}
tab <- do.call(rbind, rows)
if (is.null(tab)) stop("No file was found in either run for the groups asked: ", opt$groups)
# The LD scores of single chromosomes may exist in one run only, but only beside a
# combined table that both runs have and that was compared without a difference.
per_unit <- which(tab$kind == "per-chromosome result")
for (r in per_unit) {
  combined <- tab[tab$group == tab$group[r] & grepl("(^|/)ld_scores\\.rds$", tab$file), , drop = FALSE]
  if (nrow(combined) != 1L || grepl("^(FAIL|not compared)", combined$verdict)) {
    tab$kind[r] <- "result"
    tab$verdict[r] <- sub("^not compared, only in", "FAIL: only in", tab$verdict[r])
    tab$note[r] <- "the combined LD scores were not compared as equal in both runs, so this file counts as a result"
  }
}
# The results that a finished run of the example holds: each group has to have them
# in both runs. A pattern is met by a file that was compared (its verdict may be FAIL).
was_compared <- !grepl("^not compared", tab$verdict) & !grepl("only in [AB]$", tab$verdict)
for (group in unlist(GROUPS[asked], use.names = FALSE)) for (pattern in REQUIRED[[group]]) {
  if (any(tab$group == group & was_compared & grepl(pattern, tab$file))) next
  tab <- rbind(tab, data.frame(group = group, file = pattern, kind = "required result", size = "",
                               verdict = "FAIL: not in both runs",
                               note = "a finished run of the example holds a file of this name; no such file was compared",
                               stringsAsFactors = FALSE))
}
# The results whose number the run states itself: every chunk of the table that the
# preparation of a window scan wrote, and every model that the PI training says it
# trained. Each has to be in both runs, as the fixed names above. Without this, two
# runs that both lack the second chunk, or four of five models, would not differ.
stated_results <- function(group, folder) {
  out <- character(0)
  if (group == "window") {
    f <- file.path(folder, "shared", "01-genome_chunks.rds")
    if (file.exists(f)) out <- sprintf("results/scan_chunk_%04d\\.(csv|fst)$", as.integer(readRDS(f)$chunk_id))
  }
  if (group == "training")
    for (f in list.files(folder, pattern = "^training_metadata\\.rds$", recursive = TRUE, full.names = TRUE))
      out <- c(out, sprintf("%s/model_%d\\.rds$", basename(dirname(f)), seq_len(as.integer(readRDS(f)$n_models))))
  out
}
for (group in intersect(unlist(GROUPS[asked], use.names = FALSE), c("window", "training"))) {
  for (pattern in unique(unlist(lapply(c(fa[[group]], fb[[group]]), function(folder) stated_results(group, folder))))) {
    if (any(tab$group == group & was_compared & grepl(pattern, tab$file))) next
    tab <- rbind(tab, data.frame(group = group, file = pattern, kind = "required result", size = "",
                                 verdict = "FAIL: not in both runs",
                                 note = "the run states this result (its chunk table, or the number of models it trained); no such file was compared",
                                 stringsAsFactors = FALSE))
  }
}
# A checksum that differs is not a difference of the results when the input file
# itself was compared in this report and found to differ in run values only.
for (pc in pending_checksums) {
  explained <- vapply(pc$differing, function(input) {
    # The input was compared in this report: a row whose two files are the very bytes that the
    # two records name (their checksums are those of the records), found equal apart from run values.
    hit <- which(basename(tab$file[seq_along(paths)]) == input & tab$kind[seq_along(paths)] %in% c("R object", "GDS file"))
    hit <- hit[vapply(hit, function(r) identical(unname(tools::md5sum(paths[[r]][["a"]])), unname(pc$sums_a[[input]])) &&
                        identical(unname(tools::md5sum(paths[[r]][["b"]])), unname(pc$sums_b[[input]])), logical(1))]
    length(hit) >= 1L && all(!grepl("^FAIL", tab$verdict[hit])) && all(grepl("left out as run values", tab$note[hit]))
  }, logical(1))
  if (all(explained)) {
    tab$verdict[pc$row] <- "equal content"
    tab$note[pc$row] <- paste0(tab$note[pc$row], ": ", paste(pc$differing, collapse = ", "),
                               ", each compared above and differing in run values only")
  } else {
    tab$note[pc$row] <- paste0(tab$note[pc$row], ": ", paste(pc$differing[!explained], collapse = ", "),
                               " (not shown to differ in run values only; compare the groups training and data_prep)")
  }
}
run_tab <- if (length(run_rows)) do.call(rbind, run_rows) else NULL
if (!is.null(run_tab)) {
  about_b <- list(form = opt[["b-form"]], dir = normalizePath(opt$b), folders = fb,
                  bundle = bundle_path(opt$b, opt[["b-form"]], opt[["b-root"]]))
  run_tab$check_on_b <- mapply(check_run_value, run_tab$key, run_tab$b, run_tab$group, MoreArgs = list(b = about_b))
  run_tab$a <- substr(ifelse(is.na(run_tab$a), "", run_tab$a), 1L, 60L)
  run_tab$b <- substr(ifelse(is.na(run_tab$b), "", run_tab$b), 1L, 60L)
}

#################### THE REPORT ####################

compared <- tab[!grepl("^not compared", tab$verdict), , drop = FALSE]
n_fail <- sum(grepl("^FAIL", tab$verdict))
n_run_fail <- if (is.null(run_tab)) 0L else sum(run_tab$check_on_b == toupper(run_tab$check_on_b))   # the words in capitals
pass <- n_fail == 0L && n_run_fail == 0L && nrow(compared) > 0L
md_table <- function(d) c(paste0("| ", paste(names(d), collapse = " | "), " |"),
                          paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
                          apply(d, 1L, function(r) paste0("| ", paste(gsub("\\|", "/", r), collapse = " | "), " |")))
lines <- c(
  "# Comparison of two runs of the shipped example", "",
  sprintf("- Run A (%s form): `%s`", opt[["a-form"]], normalizePath(opt$a)),
  sprintf("- Run B (%s form): `%s`", opt[["b-form"]], normalizePath(opt$b)),
  sprintf("- Groups: %s", opt$groups),
  sprintf("- Numbers: %s", if (TOLERANCE > 0) sprintf("equal within the relative tolerance %g; everything else as text", TOLERANCE) else
            "every cell identical as text (no tolerance)"),
  sprintf("- Written by `data-example/compare-example-runs.R` on %s.", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), "",
  "## Summary", "",
  sprintf("- Files compared: %d (%d tables, %d R objects, %d settings records, %d checksum records, %d GDS files).",
          nrow(compared), sum(compared$kind == "table"), sum(compared$kind == "R object"),
          sum(compared$kind == "settings record"), sum(compared$kind == "checksum record"), sum(compared$kind == "GDS file")),
  sprintf("- Files that differ, exist in one run only, or are required and were not compared: %d.", n_fail),
  sprintf("- GLOWr version that each run recorded: A %s, B %s.", RUN_VERSION$a, RUN_VERSION$b),
  sprintf("- Records and plots counted and not compared: %d.", sum(grepl("^not compared", tab$verdict))),
  sprintf("- Values that describe the run, checked on run B: %d, of which %d fail their check.",
          if (is.null(run_tab)) 0L else nrow(run_tab), n_run_fail),
  sprintf("- Verdict: **%s**.", if (pass) "PASS" else "FAIL"), "",
  "A table is identical when it has the same columns and rows and every cell is the same text,",
  "so a missing cell is missing in both runs. An R object is identical by `identical()`. An",
  "object that holds a function is equal when `all.equal()` with tolerance 0 finds no difference.",
  "With a tolerance, a number in a table or in an R object is compared with its partner by the",
  "relative difference |a - b| / max(|a|, |b|), one pair at a time, never as a mean over a vector.",
  "A GDS file has identical content when its nodes are the same and each holds identical data.",
  "A checksum record is equal when each run's record states the checksums of the files it names",
  "and the two runs name the same inputs with the same checksums, or with checksums of files",
  "that this report compared and found to differ in run values only. The version of GLOWr in a",
  "file is a run value when it is the version that the run recorded for itself.",
  if (BUNDLE$label_only) c("", "The phenotype bundles of the two runs differ in the label of their run (`covar_set`) and in",
                           "nothing else, so the checksum that the calibration report quotes is a run value.") else NULL,
  "",
  "## Files", "", md_table(tab), "")
if (!is.null(run_tab)) lines <- c(lines,
  "## Values that describe the run", "",
  "These keys and fields name the run, its time, its software or a folder. They are left out of",
  "the comparison above and checked on run B.", "", md_table(run_tab), "")
if (!is.na(opt$out)) { writeLines(lines, opt$out); cat("Report written to", opt$out, "\n") } else cat(lines, sep = "\n")
if (n_fail > 0L) { cat("\nFailing files:\n"); print(tab[grepl("^FAIL", tab$verdict), c("group", "file", "verdict", "note")], row.names = FALSE) }
cat(sprintf("Compared %d file(s): %d differ. Run values failing their check: %d.\n", nrow(compared), n_fail, n_run_fail))
cat(sprintf("Verdict for the comparison: %s\n", if (pass) "PASS" else "FAIL"))
quit(save = "no", status = if (pass) 0L else 1L)
