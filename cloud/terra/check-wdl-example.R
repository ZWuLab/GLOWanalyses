#!/usr/bin/env Rscript
# ============================================================================
# cloud/terra/check-wdl-example.R - the WDL route's self-test
#
# Compares what glow_scan.wdl wrote for the synthetic example (run-example-wdl.sh, seven runs
# under --wdl) with what the templates wrote through run-example.sh in a reference tree
# (--reference, a GLOWanalyses root). The two routes run the same scripts in the same image,
# so the tables must be equal: the single-variant scan, its LD scores and its calibrated table,
# the per-chromosome and aggregated gene tables, the aggregated window tables, and the
# per-chromosome and aggregated coding tables. Two runs are checked for their mechanics only,
# because run-example.sh has no counterpart for them: the calibrated gene run (its Stage-1
# snapshot carries the single-variant run's inflation factor, and only GLOW's columns move)
# and the coding run with the STAAR comparator (its merged table carries STAAR_O_glowG). The
# gene run calibrated from the single-variant run's report must equal the one calibrated from
# its two files, byte for byte. The named outputs of every run are checked on disk (miniwdl
# lays each output out under out/<output name>/), with the resource monitor's exit line and
# the run_info table.
# A table is equal when its bytes are, else when its column names, dimensions and every
# value are. --tolerance T (default 0, exact) relaxes the value comparison of the two
# full-precision artifacts, ld_scores.rds and the calibrated single-variant table, for a
# reference computed in another environment (the p-value tables stay exact). Non-zero exit
# on any failure.
#
# Usage: Rscript cloud/terra/check-wdl-example.R --wdl <OUTDIR> --reference <GLOWanalyses root> [--tolerance T]

args <- commandArgs(trailingOnly = TRUE)
wdl_dir <- NA_character_; ref_dir <- NA_character_; tol <- 0
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--wdl" && i < length(args)) { wdl_dir <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--reference" && i < length(args)) { ref_dir <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--tolerance" && i < length(args)) { tol <- as.numeric(args[i + 1L]); i <- i + 2L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(wdl_dir) || is.na(ref_dir)) stop("Usage: --wdl <OUTDIR> --reference <GLOWanalyses root> [--tolerance T]")
stopifnot(is.finite(tol), tol >= 0)
stopifnot(dir.exists(wdl_dir), dir.exists(ref_dir))

.fail <- character(0)
ok <- function(cond, msg) {
  status <- isTRUE(cond)
  cat(sprintf("  [%s] %s\n", if (status) "PASS" else "FAIL", msg))
  if (!status) .fail <<- c(.fail, msg)
  invisible(status)
}
# A WDL output file, found by name under <run>/out/ (miniwdl lays outputs out by output name).
wdl_file <- function(run, name) {
  hits <- list.files(file.path(wdl_dir, run, "out"), pattern = paste0("^", gsub(".", "\\.", name, fixed = TRUE), "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0L) return(NA_character_)
  hits[1]
}
# A chr column is read as character so that "1" and "01" style values are not coerced.
read_csv <- function(f) {
  hdr <- names(read.csv(f, nrows = 1L, check.names = FALSE))
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE,
           colClasses = if ("chr" %in% hdr) c(chr = "character") else NA)
}
# Equal tables: identical bytes first (the strictest), else the same column names and
# dimensions and equal values column by column (all.equal alone ignores the names).
same_table <- function(wdl, ref, label, tolerance = 0) {
  if (is.na(wdl) || !file.exists(wdl)) return(ok(FALSE, sprintf("%s: WDL output missing", label)))
  if (!file.exists(ref)) return(ok(FALSE, sprintf("%s: reference missing (%s)", label, ref)))
  if (unname(tools::md5sum(wdl)) == unname(tools::md5sum(ref))) return(ok(TRUE, sprintf("%s: byte-identical", label)))
  a <- read_csv(wdl); b <- read_csv(ref)
  same_shape <- identical(names(a), names(b)) && identical(dim(a), dim(b))
  if (!same_shape) {
    cat(sprintf("        %d x %d vs %d x %d; column names %s\n", nrow(a), ncol(a), nrow(b), ncol(b),
                if (identical(names(a), names(b))) "equal" else "differ"))
    return(ok(FALSE, sprintf("%s: equal names, dimensions and values", label)))
  }
  eq <- isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = tolerance))
  if (!eq) {
    cols <- names(a)[vapply(names(a), function(cn) !isTRUE(all.equal(a[[cn]], b[[cn]], tolerance = tolerance)), logical(1))]
    cat(sprintf("        rows %d vs %d; differing columns: %s\n", nrow(a), nrow(b),
                if (length(cols)) paste(cols, collapse = ", ") else "(none; attributes only)"))
  }
  ok(eq, sprintf("%s: equal names, dimensions and values%s (not byte-identical)", label,
                 if (tolerance > 0) sprintf(" within %g", tolerance) else ""))
}
# A key/value table (single_variant_calibration.tsv, 01-run_info.tsv) as a named character vector.
read_kv <- function(f) { d <- read.delim(f, stringsAsFactors = FALSE, colClasses = "character"); setNames(d$value, d$key) }
# Files under a miniwdl output directory out/<output name>/ (recursive: arrays are laid out in
# numbered subdirectories); character(0) when the output was null.
out_files <- function(run, name) list.files(file.path(wdl_dir, run, "out", name), recursive = TRUE, full.names = TRUE)
factor_from_txt <- function(f) {
  ln <- grep("inflation factor", readLines(f), value = TRUE)
  as.numeric(sub(".*:\\s*([0-9.eE+-]+).*", "\\1", ln[1]))
}

cat("WDL route:", wdl_dir, "\nReference:", ref_dir, "\nTolerance for the two full-precision artifacts:", tol, "\n")

# ---- 1. The single-variant run against runs/example/single-variant ----
cat("1. single_variant\n")
sv_ref <- file.path(ref_dir, "runs/example/single-variant/outputs/results")
same_table(wdl_file("single-variant", "marginal_all.csv"), file.path(sv_ref, "marginal_all.csv"), "marginal_all.csv")
same_table(wdl_file("single-variant", "single_variant_all_calibrated.csv"),
           file.path(sv_ref, "single_variant_all_calibrated.csv"), "single_variant_all_calibrated.csv", tolerance = tol)
ld_w <- wdl_file("single-variant", "ld_scores.rds")
ld_r <- file.path(ref_dir, "data-example/cohort/example_ld/ld_scores.rds")
if (!is.na(ld_w) && file.exists(ld_r)) {
  if (tol == 0) ok(identical(readRDS(ld_w), readRDS(ld_r)), "ld_scores.rds: identical() to the reference")
  else ok(isTRUE(all.equal(readRDS(ld_w), readRDS(ld_r), tolerance = tol)), sprintf("ld_scores.rds: all.equal() to the reference within %g", tol))
} else ok(FALSE, "ld_scores.rds: WDL output or reference missing")
cal_w <- wdl_file("single-variant", "single_variant_calibration.txt")
cal_r <- file.path(sv_ref, "single_variant_calibration.txt")
f_w <- if (!is.na(cal_w)) factor_from_txt(cal_w) else NA_real_
f_r <- if (file.exists(cal_r)) factor_from_txt(cal_r) else NA_real_
ok(is.finite(f_w) && is.finite(f_r) && f_w == f_r,
   sprintf("single_variant_calibration.txt: inflation factor %s equals the reference %s", f_w, f_r))
# The calibration report: the one-file calibration input of a set-based run.
rep_w <- wdl_file("single-variant", "single_variant_calibration.tsv")
rep_factor <- NA_real_
if (!is.na(rep_w)) {
  kv <- read_kv(rep_w)
  need <- c("method", "factor", "n_variants", "pheno_bundle_md5", "chromosomes", "GLOWr", "created_at")
  rep_factor <- suppressWarnings(as.numeric(kv["factor"]))
  ok(all(need %in% names(kv)) && kv["method"] == "ldsc_intercept" && is.finite(rep_factor) && nchar(kv["pheno_bundle_md5"]) == 32L,
     sprintf("single_variant_calibration.tsv: the keys %s present; method %s; factor %s; a 32-character pheno_bundle_md5",
             paste(need, collapse = ", "), kv["method"], kv["factor"]))
  ok(is.finite(rep_factor) && is.finite(f_w) && abs(rep_factor - f_w) < 5e-5,
     sprintf("single_variant_calibration.tsv: factor %.6f equals the text file's %.4f to the printed precision", rep_factor, f_w))
} else ok(FALSE, "single_variant_calibration.tsv: WDL output missing")

# ---- 2. The gene run without calibration against runs/example/snv-set ----
cat("2. gene (no calibration)\n")
g_ref <- file.path(ref_dir, "runs/example/snv-set/outputs")
for (chr in c(21L, 22L))
  same_table(wdl_file("gene", sprintf("glow_chr%d.csv", chr)), file.path(g_ref, "results", sprintf("glow_chr%d.csv", chr)),
             sprintf("results/glow_chr%d.csv", chr))
for (f in c("glow_results_all.csv", "glow_top_hits.csv"))
  same_table(wdl_file("gene", f), file.path(g_ref, "aggregated", f), file.path("aggregated", f))

# ---- 3. The calibrated gene run: the factor travels, only GLOW's columns move ----
cat("3. gene (calibrated by the single_variant run)\n")
sh <- wdl_file("gene-calibrated", "shared.tar.gz")
if (!is.na(sh)) {
  td <- tempfile(); dir.create(td); untar(sh, exdir = td)
  cu <- readRDS(file.path(td, "run/outputs/shared/01-config_used.rds"))
  ok(is.finite(cu$inflation_factor) && is.finite(f_w) && abs(cu$inflation_factor - f_w) < 5e-5,
     sprintf("Stage-1 snapshot inflation_factor %.6f equals the single-variant factor %.4f (to the printed precision)",
             cu$inflation_factor, f_w))
  a <- wdl_file("gene-calibrated", "glow_chr21.csv"); b <- wdl_file("gene", "glow_chr21.csv")
  if (!is.na(a) && !is.na(b)) {
    ca <- read_csv(a); cb <- read_csv(b)
    glow_cols  <- grep("^GLOW_", names(ca), value = TRUE)
    other_cols <- setdiff(names(ca), c(glow_cols, "min_p_standard", "min_p_SPA"))
    ok(identical(names(ca), names(cb)) && isTRUE(all.equal(ca[other_cols], cb[other_cols], check.attributes = FALSE, tolerance = 0)),
       "glow_chr21.csv: every column other than GLOW_*, min_p_standard and min_p_SPA equals the uncalibrated run")
    ok(cu$inflation_factor == 1 || !isTRUE(all.equal(ca[glow_cols], cb[glow_cols], tolerance = 0)),
       "glow_chr21.csv: GLOW's columns differ from the uncalibrated run (the factor was applied)")
  } else ok(FALSE, "glow_chr21.csv of the calibrated or the plain gene run missing")
} else ok(FALSE, "gene-calibrated: shared.tar.gz missing")

# ---- 3b. The gene run calibrated from the report equals the one calibrated from the two files ----
cat("3b. gene (calibrated from the single_variant run's report)\n")
sh2 <- wdl_file("gene-calibrated-report", "shared.tar.gz")
if (!is.na(sh2) && !is.na(sh)) {
  td2 <- tempfile(); dir.create(td2); untar(sh2, exdir = td2)
  cu2 <- readRDS(file.path(td2, "run/outputs/shared/01-config_used.rds"))
  ok(identical(cu2$calibration$source, "report") && identical(cu2$calibration$method, "ldsc_intercept"),
     sprintf("Stage-1 snapshot: calibration source '%s', method '%s'", cu2$calibration$source, cu2$calibration$method))
  ok(is.finite(cu2$inflation_factor) && is.finite(rep_factor) && cu2$inflation_factor == rep_factor,
     sprintf("Stage-1 snapshot inflation_factor %.17g equals the report's %.17g exactly", cu2$inflation_factor, rep_factor))
  ok(is.finite(cu2$inflation_factor) && is.finite(cu$inflation_factor) && cu2$inflation_factor == cu$inflation_factor,
     sprintf("Stage-1 snapshot inflation_factor equals the two-file run's %.17g exactly", cu$inflation_factor))
  for (f in c("glow_chr21.csv", "glow_chr22.csv", "glow_results_all.csv", "glow_top_hits.csv"))
    same_table(wdl_file("gene-calibrated-report", f), wdl_file("gene-calibrated", f), paste0(f, " (report vs two files)"))
} else ok(FALSE, "gene-calibrated-report or gene-calibrated: shared.tar.gz missing")

# ---- 4. The window run against runs/example/snv-set-window ----
cat("4. window (no calibration)\n")
w_ref <- file.path(ref_dir, "runs/example/snv-set-window/outputs/aggregated")
for (f in c("scan_results_all.csv", "scan_sig_windows.csv", "scan_loci.csv"))
  same_table(wdl_file("window", f), file.path(w_ref, f), file.path("aggregated", f))

# ---- 5. The coding run against runs/example/snv-set-coding ----
cat("5. coding (no calibration, no STAAR)\n")
c_ref <- file.path(ref_dir, "runs/example/snv-set-coding/outputs")
for (chr in c(21L, 22L))
  same_table(wdl_file("coding", sprintf("glow_coding_chr%d.csv", chr)),
             file.path(c_ref, "results-glow", sprintf("glow_coding_chr%d.csv", chr)),
             sprintf("results-glow/glow_coding_chr%d.csv", chr))
for (f in c("glow_coding_all.csv", "coding_merged.csv", "coding_sig_by_category.csv", "coding_sig_by_category_cells.csv"))
  same_table(wdl_file("coding", f), file.path(c_ref, "aggregated", f), file.path("aggregated", f))

# ---- 6. The coding run with the STAAR comparator: mechanics (run-example.sh has no counterpart) ----
cat("6. coding with staar_enabled (mechanics)\n")
cs <- wdl_file("coding-staar", "coding_merged.csv")
if (!is.na(cs)) {
  d <- read_csv(cs)
  ok("STAAR_O_glowG" %in% names(d) && any(is.finite(d$STAAR_O_glowG)),
     sprintf("coding_merged.csv: %d rows; STAAR_O_glowG present with %d finite value(s)",
             nrow(d), if ("STAAR_O_glowG" %in% names(d)) sum(is.finite(d$STAAR_O_glowG)) else 0L))
  plain <- wdl_file("coding", "coding_merged.csv")
  if (!is.na(plain)) {
    e <- read_csv(plain)
    shared_cols <- intersect(names(d), names(e)); shared_cols <- shared_cols[!grepl("^(STAAR|ACAT)_", shared_cols)]
    ok(nrow(d) == nrow(e) && isTRUE(all.equal(d[shared_cols], e[shared_cols], check.attributes = FALSE, tolerance = 0)),
       sprintf("coding_merged.csv: the %d non-STAAR columns equal the run without the comparator", length(shared_cols)))
  }
} else ok(FALSE, "coding-staar: coding_merged.csv missing")

# ---- 7. The named outputs, as miniwdl laid them out (out/<output name>/), per run type ----
cat("7. named outputs on disk\n")
has_n <- function(run, name, n, what = "file(s)") {
  fs <- out_files(run, name)
  ok(length(fs) == n, sprintf("%s: output %s has %d %s (expected %d)", run, name, length(fs), what, n))
  invisible(fs)
}
has_file <- function(run, name, base) {
  fs <- out_files(run, name)
  ok(length(fs) == 1L && basename(fs[1]) == base, sprintf("%s: output %s is %s", run, name, base))
  invisible(fs)
}
absent <- function(run, name) ok(length(out_files(run, name)) == 0L, sprintf("%s: output %s is null (not this run type's)", run, name))
monitor_exit <- function(run, name, n) {
  fs <- out_files(run, name)
  good <- vapply(fs, function(f) any(grepl("^\\[monitor exit ", readLines(f, warn = FALSE))), logical(1))
  ok(length(fs) == n && all(good), sprintf("%s: output %s has %d resource file(s), each ending in the monitor's exit line", run, name, length(fs)))
}
# single_variant
has_file("single-variant", "single_variant_results", "marginal_all.csv")
has_file("single-variant", "single_variant_calibrated", "single_variant_all_calibrated.csv")
has_file("single-variant", "calibration_report_out", "single_variant_calibration.tsv")
has_file("single-variant", "single_variant_calibration", "single_variant_calibration.txt")
has_file("single-variant", "ld_scores_out", "ld_scores.rds")
has_file("single-variant", "single_variant_plots", "marginal_plots.pdf")
has_file("single-variant", "calibrate_log", "calibrate.log")
has_n("single-variant", "per_chrom_logs", 2L)
monitor_exit("single-variant", "per_chrom_resources", 2L)
monitor_exit("single-variant", "calibrate_resources", 1L)
absent("single-variant", "combined_results"); absent("single-variant", "run_info")
# gene
has_file("gene", "combined_results", "glow_results_all.csv")
has_file("gene", "top_hits", "glow_top_hits.csv")
ri <- has_file("gene", "run_info", "01-run_info.tsv")
if (length(ri) == 1L) {
  kv <- read_kv(ri)
  ok(identical(unname(kv["region_type"]), "gene") && identical(unname(kv["calibration_source"]), "none") &&
       identical(unname(kv["chromosomes"]), "21,22") && nzchar(kv["GLOWr"]),
     sprintf("gene: run_info says region_type %s, calibration_source %s, chromosomes %s, GLOWr %s",
             kv["region_type"], kv["calibration_source"], kv["chromosomes"], kv["GLOWr"]))
}
ri2 <- out_files("gene-calibrated-report", "run_info")
if (length(ri2) == 1L) { kv2 <- read_kv(ri2); ok(identical(unname(kv2["calibration_source"]), "report"), sprintf("gene-calibrated-report: run_info says calibration_source %s", kv2["calibration_source"])) }
has_n("gene", "per_chrom_results", 2L)
has_n("gene", "per_chrom_logs", 2L)
monitor_exit("gene", "per_chrom_resources", 2L)
monitor_exit("gene", "prepare_resources", 1L)
monitor_exit("gene", "aggregate_resources", 1L)
ok(length(out_files("gene", "plots")) >= 1L && all(grepl("\\.pdf$", out_files("gene", "plots"))),
   sprintf("gene: output plots has %d pdf(s)", length(out_files("gene", "plots"))))
has_file("gene", "aggregated_tar", "aggregated.tar.gz"); has_file("gene", "shared_artifacts", "shared.tar.gz")
has_file("gene", "prepare_log", "prepare.log"); has_file("gene", "aggregate_log", "aggregate.log")
absent("gene", "sig_windows"); absent("gene", "coding_merged"); absent("gene", "noncoding_completeness"); absent("gene", "single_variant_results")
# window
has_file("window", "combined_results", "scan_results_all.csv")
has_file("window", "sig_windows", "scan_sig_windows.csv"); has_file("window", "loci", "scan_loci.csv")
absent("window", "top_hits")
# coding
has_file("coding", "combined_results", "glow_coding_all.csv")
has_file("coding", "coding_merged", "coding_merged.csv"); has_file("coding", "coding_sig_by_category", "coding_sig_by_category.csv")
has_file("coding", "coding_sig_cells", "coding_sig_by_category_cells.csv"); has_file("coding", "lambda_by_category", "lambda_by_category.csv")
has_n("coding", "per_chrom_results", 2L)
absent("coding", "top_hits"); absent("coding", "noncoding_merged")

cat("\n")
if (length(.fail)) {
  cat(sprintf("%d CHECK(S) FAILED:\n", length(.fail))); cat(paste0("  - ", .fail, "\n"), sep = "")
  quit(save = "no", status = 1L)
}
cat("All checks passed: the WDL route equals the template route on the synthetic example.\n")
