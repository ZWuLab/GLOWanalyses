#!/usr/bin/env Rscript
# ============================================================================
# 02-single-variant/calibrate.R - genomic-control calibration
#
# Cohort-agnostic genomic-control calibration step. COMPUTES this run's inflation
# factor from its scan (results/marginal_all.csv)
# + the shared LD scores (ld_scores_path) via GLOWr::estimate_inflation_factor(),
# then rescales the scan's chi-squares by it (chi2_cal = chi2 / factor) on the
# P-VALUE columns (pvalue [+ pvalue_SPA]); the Z / score columns are kept as the
# raw single-variant statistics (the calibration is on the chi-square scale).
#
#   method "ldsc_intercept" (default): LD Score regression intercept - the
#     confounding-only factor (needs the LD scores); preserves polygenic signal.
#   method "lambda_gc": the genomic-control lambda (no LD scores).
#
# Inputs:
#   <output_dir>/results/marginal_all.csv   (this run's scan; marginal-scan.R)
#   ld_scores_path                          (shared; compute-ld-scores.R; ldsc only)
# Outputs (in this run's results/):
#   single_variant_all_calibrated.csv   (marginal_all.csv schema; calibrated p-values)
#   single_variant_calibration.txt      (provenance: factor, method, raw->calibrated lambda)
#   single_variant_calibration.tsv      (the same facts as a key/value table, the factor at
#                                        full precision: the one-file calibration input of a
#                                        set-based run, 03-snv-set/prepare.R calibration$report)
#
# Usage (from project root; run compute-ld-scores.R first for ldsc_intercept):
#   conda activate r_env
#   Rscript 02-single-variant/calibrate.R \
#       --config <run>/config.R [--method ldsc_intercept|lambda_gc]

suppressMessages({
  library(data.table)
  library(GLOWr)   # estimate_inflation_factor,
})                                                       # calibrate_pvalues, compute_lambda_gc
null_or <- function(a, b) if (is.null(a)) b else a

# ---- CSV writer: data.table::fwrite() unless a subnormal double is present ----
# The helper of marginal-scan.R: fwrite() mis-writes doubles below .Machine$double.xmin
# (1e-308 came out as 1.6e-308; data.table 1.15.4 and 1.18.4; review round
# is written by write.csv(), exact and slower; every other table takes the fast writer.
write_csv_exact <- function(d, path) {
  has_subnormal <- function(x) is.double(x) &&
    any(x != 0 & abs(x) < .Machine$double.xmin, na.rm = TRUE)
  if (any(vapply(d, has_subnormal, logical(1)))) {
    cat(sprintf("  %s holds a value below %.3g: writing with write.csv (exact, slower)\n",
                basename(path), .Machine$double.xmin))
    utils::write.csv(d, path, row.names = FALSE, na = "NA")
  } else {
    data.table::fwrite(d, path, na = "NA")
  }
}

# ---- Args: --config <path> [--method M] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; method_arg <- NULL
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--method" && i < length(args)) { method_arg <- args[i + 1L]; i <- i + 2L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))

# ---- Config sourcing (the proven idiom) ----
source(config_path, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

# ---- Run identity + paths ----
output_dir  <- file.path(dirname(config_path), "outputs")
run_name    <- basename(dirname(config_path))
# Shared-output OPT-IN: NULL (default) keeps outputs/ local; a shared-root path
# born-shares this run's outputs/ (symlink into <root>/<repo-relative-path>) before it
# is populated. Normally marginal-scan.R already created/born-shared it (then this is a
# no-op); the guard makes calibrate self-sufficient if run first. GLOWpipeline is loaded
ssr <- get0("symlinked_shared_root", envir = cfg, ifnotfound = NULL, inherits = FALSE)
if (!is.null(ssr)) {
  suppressMessages(library(GLOWpipeline))
  GLOWpipeline::ensure_shared_output_dir(output_dir, share_root = ssr)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
results_dir <- file.path(output_dir, "results")
marg_path   <- file.path(results_dir, "marginal_all.csv")        # this run's scan
ld_path     <- g0("ld_scores_path")                              # shared (compute-ld-scores.R)
out_csv     <- file.path(results_dir, "single_variant_all_calibrated.csv")
out_prov    <- file.path(results_dir, "single_variant_calibration.txt")
method      <- null_or(method_arg, null_or(g0("calibration_method"), "ldsc_intercept"))

if (!file.exists(marg_path))
  stop("Scan not found: ", marg_path, "\nRun marginal-scan.R (default / --combine) first.")

# ---- 1. Compute the inflation factor for this run (self-sustained) ----
d <- data.table::fread(marg_path, colClasses = list(character = "chr"), data.table = FALSE)   # one copy: as.data.frame() on a data.table copies it again
ld_scores <- NULL
if (method == "ldsc_intercept") {
  if (is.null(ld_path) || !file.exists(ld_path)) {
    stop("method 'ldsc_intercept' needs LD scores at ", null_or(ld_path, "<unset ld_scores_path>"),
         "\nRun compute-ld-scores.R first (or use --method lambda_gc).")
  }
  ld_scores <- readRDS(ld_path)
}
est    <- estimate_inflation_factor(d, method = method, ld_scores = ld_scores)
factor <- est$factor
if (!is.finite(factor) || factor <= 0) {
  stop(sprintf("estimated inflation factor is not usable (%.4g); aborting.", factor))
}
cat(sprintf("[%s] M1 '%s' (%s): inflation factor = %.4f (n=%d)\n",
            format(Sys.time(), "%H:%M:%S"), run_name, method, factor, est$n_variants))

# ---- 2. Calibrate the single-variant p-values: chi2 -> chi2 / factor ----
# The headline lambdas are computed on variants with MAC >= lambda_mac_min (config; default
# 20). On sequencing data most variants are ultra-rare and their p-values, standard or SPA,
# take a few discrete values (a singleton has one of two), so a median over them says nothing
# about calibration in either direction. The per-MAC-range medians of the raw p-values are
# tabulated beside the headline, so a reader sees where a genome-wide lambda comes from.
lambda_mac_min <- as.integer(null_or(g0("lambda_mac_min"), 20L))
stopifnot(length(lambda_mac_min) == 1L, is.finite(lambda_mac_min), lambda_mac_min >= 0L)
if (!"MAC" %in% names(d)) stop("marginal_all.csv has no MAC column; the lambda report needs it.")
has_spa  <- "pvalue_SPA" %in% names(d)
head_sel <- is.finite(d$MAC) & d$MAC >= lambda_mac_min
lambda_or_na <- function(p) { p <- p[is.finite(p) & p > 0 & p <= 1]; if (length(p)) compute_lambda_gc(p) else NA_real_ }
mac_labels <- c("1", "2-5", "6-19", "20-99", "100+")
mac_bin    <- cut(d$MAC, breaks = c(1, 2, 6, 20, 100, Inf), right = FALSE, labels = mac_labels)
lambda_by_mac <- data.frame(
  mac_range       = mac_labels,
  n               = as.integer(table(factor(mac_bin, levels = mac_labels))),
  lambda_standard = vapply(mac_labels, function(b) lambda_or_na(d$pvalue[!is.na(mac_bin) & mac_bin == b]), numeric(1)),
  lambda_SPA      = if (has_spa) vapply(mac_labels, function(b) lambda_or_na(d$pvalue_SPA[!is.na(mac_bin) & mac_bin == b]), numeric(1)) else NA_real_,
  stringsAsFactors = FALSE, row.names = NULL)
raw_lambda     <- lambda_or_na(d$pvalue[head_sel])
d$pvalue       <- calibrate_pvalues(d$pvalue, method = method,
                                    calibration_factor = factor)$p
cal_lambda     <- lambda_or_na(d$pvalue[head_sel])
raw_lambda_spa <- cal_lambda_spa <- NA_real_
if (has_spa) {
  raw_lambda_spa <- lambda_or_na(d$pvalue_SPA[head_sel])
  d$pvalue_SPA   <- calibrate_pvalues(d$pvalue_SPA, method = method,
                                      calibration_factor = factor)$p
  cal_lambda_spa <- lambda_or_na(d$pvalue_SPA[head_sel])
}
write_csv_exact(d, out_csv)

# ---- 3. Provenance + data-derived summary (descriptive labels; a median ratio establishes
# neither that confounding was removed nor that the tails are calibrated) ----
verdict <- function(lam) {
  if (is.na(lam)) "NA"
  else if (abs(lam - 1) <= 0.05) sprintf("%.4f (within 0.05 of 1)", lam)
  else if (lam > 1) sprintf("%.4f (above 1)", lam)
  else sprintf("%.4f (below 1)", lam)
}
fmt_lam <- function(x) ifelse(is.na(x), "NA", sprintf("%.4f", x))
mac_table_lines <- c(
  sprintf("  %-10s %12s %16s %16s", "MAC range", "variants", "lambda standard", "lambda SPA"),
  sprintf("  %-10s %12d %16s %16s", lambda_by_mac$mac_range, lambda_by_mac$n,
          fmt_lam(lambda_by_mac$lambda_standard), fmt_lam(lambda_by_mac$lambda_SPA)),
  "  (NA: fewer than 50 usable p-values in the range, finite and in (0, 1], the minimum compute_lambda_gc() accepts)")
git_sha <- tryCatch(suppressWarnings(system("git rev-parse --short HEAD",
                    intern = TRUE, ignore.stderr = TRUE)), error = function(e) NA)
if (length(git_sha) != 1L || !nzchar(git_sha)) git_sha <- "unknown"

prov <- c(
  sprintf("Single-variant (M1) calibration - config '%s' (%s)",
          run_name, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("  method            : %s", method),
  sprintf("  inflation factor  : %.4f   (estimate_inflation_factor, n=%d)",
          factor, est$n_variants),
  if (method == "ldsc_intercept") c(
    sprintf("  LDSC intercept    : %.4f (SE %.4f); slope %.4g (SE %.4g); mean chi-square %.4f",
            est$intercept, est$intercept_se, est$slope, est$slope_se, est$mean_chi2),
    sprintf("  confounding ratio : %.3f (SE %.3f)   ((intercept - 1) / (mean chi-square - 1); LD scores from compute-ld-scores.R)",
            est$confounding_ratio, est$ratio_se))
  else "  confounding ratio : (lambda_gc removes all inflation, not just confounding)",
  sprintf("  git SHA           : %s", git_sha),
  sprintf("  variants (scan)   : %d", nrow(d)),
  "",
  sprintf("lambda_GC of the single-variant scan on the %d variants with MAC >= %d (raw -> calibrated):",
          sum(head_sel), lambda_mac_min),
  sprintf("  pvalue (standard) : %s  ->  %s", verdict(raw_lambda),     verdict(cal_lambda)),
  if (has_spa)
  sprintf("  pvalue_SPA        : %s  ->  %s", verdict(raw_lambda_spa), verdict(cal_lambda_spa))
  else "  pvalue_SPA        : (not present)",
  "",
  "lambda_GC of the raw p-values per MAC range (every scanned variant):",
  mac_table_lines,
  "",
  paste0("NOTE: lambda_GC is the median chi-square over the variants it is computed on, divided by ",
         "the null median. It summarizes the center of the p-value distribution, not its tails, ",
         "and on ultra-rare variants (MAC in the single digits) the p-values are discrete: a ",
         "singleton carried by a case or by a control has one of two values, standard or SPA ",
         "(an intercept-only binary null model gives them in closed form), so a median over ",
         "them is not a calibration measure in either direction. ",
         if (method == "ldsc_intercept")
           paste0("The factor is the LDSC intercept: the regression of the standard chi-squares ",
                  "on the LD scores over the variants with MAF >= 0.01 outside the MHC that have ",
                  "an LD score; it removes the part of the inflation that behaves like ",
                  "confounding and keeps polygenic signal. ")
         else
           paste0("The factor is lambda_GC itself, the median chi-square ratio of the standard ",
                  "p-values over the variants with MAF >= 0.01; it removes every inflation of ",
                  "the center, confounding or polygenic. "),
         "Calibration is on the chi-square/p-value scale, and the Z/score columns stay the raw ",
         "statistics.")
)
writeLines(prov, out_prov)
cat(paste(prov, collapse = "\n"), "\n")

# ---- 4. The calibration report: single_variant_calibration.tsv ----
# The same facts as the provenance text, as one key/value table with the factor at full
# precision, so that a set-based run can take its calibration from this one file
# (03-snv-set/prepare.R, calibration$report) instead of re-estimating it from
# marginal_all.csv and the LD scores. pheno_bundle_md5 is the checksum of the phenotype
# bundle this scan ran on: the factor is specific to its samples, outcome and covariates,
# and prepare.R refuses the report for a run on another bundle.
pheno_path <- g0("pheno_path")
pkg_ver    <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
fmt        <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) NA_character_ else sprintf("%.10g", x)
release_source <- if (file.exists(".release-source")) readLines(".release-source", n = 1L, warn = FALSE) else "development tree"
report <- c(
  method                = method,
  factor                = sprintf("%.17g", factor),            # full precision: a run calibrated from the report equals one calibrated from the two files
  n_variants            = as.character(est$n_variants),       # variants in the regression (ldsc) or the lambda (lambda_gc)
  confounding_ratio     = if (method == "ldsc_intercept") fmt(est$confounding_ratio) else NA_character_,
  ratio_se              = if (method == "ldsc_intercept") fmt(est$ratio_se) else NA_character_,
  intercept_se          = if (method == "ldsc_intercept") fmt(est$intercept_se) else NA_character_,
  slope                 = if (method == "ldsc_intercept") fmt(est$slope) else NA_character_,
  slope_se              = if (method == "ldsc_intercept") fmt(est$slope_se) else NA_character_,
  mean_chisq            = if (method == "ldsc_intercept") fmt(est$mean_chi2) else NA_character_,
  lambda_mac_min        = as.character(lambda_mac_min),       # the headline lambdas below are computed on MAC >= this
  n_variants_lambda     = as.character(sum(head_sel)),
  lambda_raw            = fmt(raw_lambda),
  lambda_calibrated     = fmt(cal_lambda),
  lambda_raw_SPA        = fmt(raw_lambda_spa),
  lambda_calibrated_SPA = fmt(cal_lambda_spa),
  setNames(as.character(lambda_by_mac$n), paste0("n_mac_", lambda_by_mac$mac_range)),
  setNames(vapply(lambda_by_mac$lambda_standard, fmt, character(1)), paste0("lambda_standard_mac_", lambda_by_mac$mac_range)),
  setNames(vapply(lambda_by_mac$lambda_SPA, fmt, character(1)), paste0("lambda_SPA_mac_", lambda_by_mac$mac_range)),
  n_variants_scan       = as.character(nrow(d)),
  chromosomes           = paste(unique(d$chr), collapse = ","),
  pheno_bundle          = if (!is.null(pheno_path)) basename(pheno_path) else NA_character_,
  pheno_bundle_md5      = if (!is.null(pheno_path) && file.exists(pheno_path)) unname(tools::md5sum(pheno_path)) else NA_character_,
  ld_scores             = if (method == "ldsc_intercept") basename(ld_path) else NA_character_,
  run_name              = run_name,
  GFisher               = pkg_ver("GFisher"),
  GLOWr                 = pkg_ver("GLOWr"),
  GLOWpipeline          = pkg_ver("GLOWpipeline"),
  GLOWanalyses          = release_source,
  created_at            = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
out_tsv <- file.path(results_dir, "single_variant_calibration.tsv")
utils::write.table(data.frame(key = names(report), value = unname(report), stringsAsFactors = FALSE),
                   out_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("\n[%s] Wrote %s, %s and %s\n", format(Sys.time(), "%H:%M:%S"), out_csv, out_prov, out_tsv))
