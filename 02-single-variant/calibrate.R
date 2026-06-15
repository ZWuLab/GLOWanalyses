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
d <- as.data.frame(data.table::fread(marg_path, colClasses = list(character = "chr")))
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
has_spa        <- "pvalue_SPA" %in% names(d)
raw_lambda     <- compute_lambda_gc(d$pvalue)
d$pvalue       <- calibrate_pvalues(d$pvalue, method = method,
                                    calibration_factor = factor)$p
cal_lambda     <- compute_lambda_gc(d$pvalue)
raw_lambda_spa <- cal_lambda_spa <- NA_real_
if (has_spa) {
  raw_lambda_spa <- compute_lambda_gc(d$pvalue_SPA)
  d$pvalue_SPA   <- calibrate_pvalues(d$pvalue_SPA, method = method,
                                      calibration_factor = factor)$p
  cal_lambda_spa <- compute_lambda_gc(d$pvalue_SPA)
}
data.table::fwrite(d, out_csv)

# ---- 3. Provenance + data-derived summary (no hardcoded conclusions) ----
verdict <- function(lam) {
  if (is.na(lam)) "NA"
  else if (abs(lam - 1) <= 0.05) sprintf("%.4f (~1: confounding removed)", lam)
  else if (lam > 1) sprintf("%.4f (residual inflation remains)", lam)
  else sprintf("%.4f (over-corrected / deflated)", lam)
}
git_sha <- tryCatch(suppressWarnings(system("git rev-parse --short HEAD",
                    intern = TRUE, ignore.stderr = TRUE)), error = function(e) NA)
if (length(git_sha) != 1L || !nzchar(git_sha)) git_sha <- "unknown"

prov <- c(
  sprintf("Single-variant (M1) calibration - config '%s' (%s)",
          run_name, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("  method            : %s", method),
  sprintf("  inflation factor  : %.4f   (estimate_inflation_factor, n=%d)",
          factor, est$n_variants),
  if (method == "ldsc_intercept")
    sprintf("  confounding ratio : %.3f   (LDSC; LD scores from compute-ld-scores.R)", est$confounding_ratio)
  else "  confounding ratio : (lambda_gc removes all inflation, not just confounding)",
  sprintf("  git SHA           : %s", git_sha),
  sprintf("  variants (scan)   : %d", nrow(d)),
  "",
  "lambda_GC of the single-variant scan (raw -> calibrated; calibrated should be ~1):",
  sprintf("  pvalue (standard) : %s  ->  %s", verdict(raw_lambda),     verdict(cal_lambda)),
  if (has_spa)
  sprintf("  pvalue_SPA        : %s  ->  %s", verdict(raw_lambda_spa), verdict(cal_lambda_spa))
  else "  pvalue_SPA        : (not present)",
  "",
  paste0("NOTE: ldsc_intercept removes the CONFOUNDING part of the inflation and ",
         "preserves polygenic signal, so a calibrated lambda slightly above 1 is ",
         "expected when the residual is partly polygenic (design M1). Calibration is on ",
         "the chi-square/p-value scale; Z/score columns are the raw statistics.")
)
writeLines(prov, out_prov)
cat(paste(prov, collapse = "\n"), "\n")
cat(sprintf("\n[%s] Wrote %s and %s\n", format(Sys.time(), "%H:%M:%S"), out_csv, out_prov))
