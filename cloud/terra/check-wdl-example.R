#!/usr/bin/env Rscript
# ============================================================================
# cloud/terra/check-wdl-example.R - the WDL route's self-test
#
# Compares what glow_scan.wdl wrote for the synthetic example (run-example-wdl.sh, six runs
# under --wdl) with what the templates wrote through run-example.sh in a reference tree
# (--reference, a GLOWanalyses root). The two routes run the same scripts in the same image,
# so the tables must be equal: the single-variant scan, its LD scores and its calibrated table,
# the per-chromosome and aggregated gene tables, the aggregated window tables, and the
# per-chromosome and aggregated coding tables. Two runs are checked for their mechanics only,
# because run-example.sh has no counterpart for them: the calibrated gene run (its Stage-1
# snapshot carries the single-variant run's inflation factor, and only GLOW's columns move)
# and the coding run with the STAAR comparator (its merged table carries STAAR_O_glowG).
# A table is equal when its bytes are, else when its column names, dimensions and every
# value are. Non-zero exit on any failure.
#
# Usage: Rscript cloud/terra/check-wdl-example.R --wdl <OUTDIR> --reference <GLOWanalyses root>

args <- commandArgs(trailingOnly = TRUE)
wdl_dir <- NA_character_; ref_dir <- NA_character_
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--wdl" && i < length(args)) { wdl_dir <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--reference" && i < length(args)) { ref_dir <- args[i + 1L]; i <- i + 2L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(wdl_dir) || is.na(ref_dir)) stop("Usage: --wdl <OUTDIR> --reference <GLOWanalyses root>")
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
same_table <- function(wdl, ref, label) {
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
  eq <- isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 0))
  if (!eq) {
    cols <- names(a)[vapply(names(a), function(cn) !isTRUE(all.equal(a[[cn]], b[[cn]], tolerance = 0)), logical(1))]
    cat(sprintf("        rows %d vs %d; differing columns: %s\n", nrow(a), nrow(b),
                if (length(cols)) paste(cols, collapse = ", ") else "(none; attributes only)"))
  }
  ok(eq, sprintf("%s: equal names, dimensions and values (not byte-identical)", label))
}
factor_from_txt <- function(f) {
  ln <- grep("inflation factor", readLines(f), value = TRUE)
  as.numeric(sub(".*:\\s*([0-9.eE+-]+).*", "\\1", ln[1]))
}

cat("WDL route:", wdl_dir, "\nReference:", ref_dir, "\n")

# ---- 1. The single-variant run against runs/example/single-variant ----
cat("1. single_variant\n")
sv_ref <- file.path(ref_dir, "runs/example/single-variant/outputs/results")
same_table(wdl_file("single-variant", "marginal_all.csv"), file.path(sv_ref, "marginal_all.csv"), "marginal_all.csv")
same_table(wdl_file("single-variant", "single_variant_all_calibrated.csv"),
           file.path(sv_ref, "single_variant_all_calibrated.csv"), "single_variant_all_calibrated.csv")
ld_w <- wdl_file("single-variant", "ld_scores.rds")
ld_r <- file.path(ref_dir, "data-example/cohort/example_ld/ld_scores.rds")
if (!is.na(ld_w) && file.exists(ld_r)) {
  ok(identical(readRDS(ld_w), readRDS(ld_r)), "ld_scores.rds: identical() to the reference")
} else ok(FALSE, "ld_scores.rds: WDL output or reference missing")
cal_w <- wdl_file("single-variant", "single_variant_calibration.txt")
cal_r <- file.path(sv_ref, "single_variant_calibration.txt")
f_w <- if (!is.na(cal_w)) factor_from_txt(cal_w) else NA_real_
f_r <- if (file.exists(cal_r)) factor_from_txt(cal_r) else NA_real_
ok(is.finite(f_w) && is.finite(f_r) && f_w == f_r,
   sprintf("single_variant_calibration.txt: inflation factor %s equals the reference %s", f_w, f_r))

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

cat("\n")
if (length(.fail)) {
  cat(sprintf("%d CHECK(S) FAILED:\n", length(.fail))); cat(paste0("  - ", .fail, "\n"), sep = "")
  quit(save = "no", status = 1L)
}
cat("All checks passed: the WDL route equals the template route on the synthetic example.\n")
