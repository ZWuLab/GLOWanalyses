#!/usr/bin/env Rscript
# ============================================================================
# 02-single-variant/marginal-scan.R - per-variant association scan
#
# Cohort-agnostic per-variant (single-variant) association scan. Reads the
# assemble_pheno_covar bundle (pheno_path), fits the null model once
# (fit_null_model), then scans each chromosome's (a)GDS with GLOWr::marginal_scan
# (per-chr CSV) and, in the combining modes, rbinds them into marginal_all.csv
# and draws the Manhattan/QQ plots. output_dir is derived from --config.
#
# Modes:
#   --chr N      scan ONLY chr N -> results/marginal_chr<N>.csv (array-task unit;
#                no combine, no plots).
#   (default)    loop `chroms`, per-chr CSV, then rbind -> results/marginal_all.csv
#                + plots (the all-chromosome path).
#   --combine    skip scanning; rbind existing results/marginal_chr*.csv (numeric
#                chr order) -> results/marginal_all.csv + plots (post-array).
#
# Output (under <output_dir> = <dir of --config>/outputs/):
#   results/marginal_chr<N>.csv   per-chromosome scan
#   results/marginal_all.csv      combined (default / --combine)
#   plots/marginal_plots.pdf      Manhattan + QQ (default / --combine)
#   logs/scan_log.txt             timestamped run log
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 02-single-variant/marginal-scan.R \
#       --config <run>/config.R [--chr 22 | --combine]

suppressMessages(library(GLOWr))
library(SeqArray)
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <path> [--chr N | --combine] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; chr <- NULL; combine_only <- FALSE
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--chr" && i < length(args)) { chr <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--combine") { combine_only <- TRUE; i <- i + 1L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))
if (!is.null(chr) && combine_only) stop("--chr and --combine are mutually exclusive.")
if (!is.null(chr)) stopifnot(chr >= 1L, chr <= 22L)

# ---- Config sourcing (the proven 03-snv-set / 00-data-prep idiom) ----
source(config_path, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

# ---- Run identity + output dirs (derived from the config path) ----
output_dir  <- file.path(dirname(config_path), "outputs")
run_name    <- basename(dirname(config_path))
results_dir <- file.path(output_dir, "results")
plots_dir   <- file.path(output_dir, "plots")
logs_dir    <- file.path(output_dir, "logs")
# Shared-output OPT-IN: NULL (default) keeps outputs/ local; a shared-root path
# born-shares this run's outputs/ (symlink into <root>/<repo-relative-path>) BEFORE
# the results/plots/logs subdirs below populate it. GLOWpipeline is loaded only on
ssr <- get0("symlinked_shared_root", envir = cfg, ifnotfound = NULL, inherits = FALSE)
if (!is.null(ssr)) {
  suppressMessages(library(GLOWpipeline))
  GLOWpipeline::ensure_shared_output_dir(output_dir, share_root = ssr)
}
for (d in c(results_dir, plots_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Required cohort fields + scan knobs (knobs from base) ----
pheno_path  <- g0("pheno_path"); gds_dir <- g0("gds_dir"); gds_pattern <- g0("gds_pattern")
if (is.null(pheno_path) || is.null(gds_dir) || is.null(gds_pattern))
  stop("Config must set `pheno_path`, `gds_dir`, and `gds_pattern`.")
chroms             <- null_or(g0("chroms"), 1:22)
use_SPA            <- isTRUE(null_or(g0("use_SPA"), TRUE))
chunk_size         <- as.integer(null_or(g0("chunk_size"), 2000L))
mac_cutoff         <- as.integer(null_or(g0("mac_cutoff"), 1L))
missing_imputation <- null_or(g0("missing_imputation"), "mean")
n_cores            <- as.integer(null_or(g0("n_cores"), 1L))   # forked workers per chromosome

# ---- Logging ----
log_con <- file(file.path(logs_dir, "scan_log.txt"), open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE, after = FALSE)
log_msg <- function(msg) {
  ts_msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
  cat(ts_msg, "\n"); writeLines(ts_msg, log_con); flush(log_con)
}

# ---- CSV writer: data.table::fwrite() unless a subnormal double is present ----
# fwrite() mis-writes doubles below .Machine$double.xmin (about 2.2e-308): 1e-308 came out as
# 1.6e-308 and the smallest positive double as 1.1e-308 (data.table 1.15.4 and 1.18.4; review
# small needs |Z| above about 37.5, rare but possible on a genome, so a table holding one is
# written by write.csv(), exact at 15 significant digits and slower; every other table takes
# the fast writer. The parsed values are the same either way. calibrate.R has the same helper.
write_csv_exact <- function(d, path) {
  has_subnormal <- function(x) is.double(x) &&
    any(x != 0 & abs(x) < .Machine$double.xmin, na.rm = TRUE)
  if (any(vapply(d, has_subnormal, logical(1)))) {
    log_msg(sprintf("  %s holds a value below %.3g: writing with write.csv (exact, slower)",
                    basename(path), .Machine$double.xmin))
    utils::write.csv(d, path, row.names = FALSE, na = "NA")
  } else {
    data.table::fwrite(d, path, na = "NA")
  }
}

# Per-chr GDS path from the cohort pattern (literal "{chr}" token).
gds_path_for <- function(k) file.path(gds_dir, gsub("{chr}", k, gds_pattern, fixed = TRUE))

t_total_start <- proc.time()
mode_label <- if (!is.null(chr)) sprintf("single-chr (chr%d)", chr) else if (combine_only) "combine-only" else "all-chr"
log_msg("=== Marginal Scan Started ===")
log_msg(sprintf("Run: %s | mode: %s", run_name, mode_label))
log_msg(sprintf("use_SPA: %s | chunk_size: %d | mac_cutoff: %d | imputation: %s | n_cores: %d",
                use_SPA, chunk_size, mac_cutoff, missing_imputation, n_cores))

# ===========================================================================
# Combine helper: bind the per-chr results (numeric chr order) -> marginal_all.csv + plots.
# data.table's fread, rbindlist and fwrite (data.table is a GLOWr dependency), because a
# whole genome is about 125 million rows: read.csv, do.call(rbind) and write.csv took nine
# minutes for 11 million rows and hold several copies at once. The per-chromosome tables
# are released before the plots are drawn.
#
# `in_memory` (default mode) is the list of marginal_scan() return data.frames.
# `--combine` passes NULL, so we re-read the per-chr CSVs that marginal_scan()
# wrote (post-array; same chr/data, the float round-trip is harmless - calibrate
# consumes values, not bytes). Every chromosome of `chroms` must have its table: a
# missing one stops the combine, so an incomplete scatter cannot produce a genome with
# a chromosome silently absent.
# ===========================================================================
combine_and_plot <- function(in_memory = NULL) {
  log_msg("Combining per-chromosome results")
  if (is.null(in_memory)) {
    chr_files <- file.path(results_dir, sprintf("marginal_chr%d.csv", chroms))
    absent    <- chr_files[!file.exists(chr_files)]
    if (length(absent) > 0L)
      stop("Per-chromosome table(s) missing for requested chromosome(s): ",
           paste(basename(absent), collapse = ", "),
           ". Every chromosome of `chroms` must be scanned before --combine.")
    results_list <- lapply(chr_files, function(f)
      data.table::fread(f, colClasses = list(character = "chr")))
  } else {
    results_list <- in_memory[!vapply(in_memory, is.null, logical(1))]
  }
  all_results <- data.table::rbindlist(results_list, use.names = TRUE)
  rm(results_list)                      # the per-chromosome copies go before the plots
  data.table::setDF(all_results)
  invisible(gc())

  csv_all <- file.path(results_dir, "marginal_all.csv")
  write_csv_exact(all_results, csv_all)
  log_msg(sprintf("  Combined: %d total variants -> %s", nrow(all_results), csv_all))

  # ---- Summary statistics ----
  log_msg(sprintf("  Variants with pvalue < 5e-8: %d",
                  sum(all_results$pvalue < 5e-8, na.rm = TRUE)))
  log_msg(sprintf("  Variants with pvalue < 1e-5: %d",
                  sum(all_results$pvalue < 1e-5, na.rm = TRUE)))

  # ---- Plots ----
  log_msg("Generating plots")
  pdf_path <- file.path(plots_dir, "marginal_plots.pdf")
  pdf(pdf_path, width = 14, height = 6)

  plot_manhattan(all_results, pvalue_col = "pvalue",
                 title = paste("Manhattan (standard) -", run_name))
  lambda_std <- plot_qq(all_results, pvalue_col = "pvalue",
                        title = paste("QQ (standard) -", run_name))

  if ("pvalue_SPA" %in% names(all_results)) {
    par(mfrow = c(1L, 2L))
    plot_qq(all_results, pvalue_col = "pvalue",
            title = paste("QQ Standard -", run_name))
    lambda_spa <- plot_qq(all_results, pvalue_col = "pvalue_SPA",
                          title = paste("QQ SPA -", run_name))
    par(mfrow = c(1L, 1L))
    plot_manhattan(all_results, pvalue_col = "pvalue_SPA",
                   title = paste("Manhattan (SPA) -", run_name))
    log_msg(sprintf("  Lambda_GC (standard): %.4f", lambda_std))
    log_msg(sprintf("  Lambda_GC (SPA):      %.4f", lambda_spa))
  } else {
    log_msg(sprintf("  Lambda_GC (standard): %.4f", lambda_std))
  }
  dev.off()
  log_msg(sprintf("  Plots saved to: %s", pdf_path))
  invisible(all_results)
}

# ===========================================================================
# --combine: skip scanning entirely (post-array path)
# ===========================================================================
if (combine_only) {
  combine_and_plot()
  log_msg(sprintf("=== Marginal Scan (combine) Complete in %.1f s ===",
                  (proc.time() - t_total_start)[3]))
  quit(save = "no", status = 0L)
}

# ===========================================================================
# Step 0: every requested chromosome's GDS must exist before any work. A missing file was
# skipped with a warning, so the default route finished with exit status 0 and a genome
# ===========================================================================
gds_paths <- vapply(if (!is.null(chr)) chr else chroms, gds_path_for, character(1))
absent    <- gds_paths[!file.exists(gds_paths)]
if (length(absent) > 0L)
  stop("GDS file(s) missing for requested chromosome(s): ", paste(absent, collapse = ", "),
       ". Every chromosome of `chroms` must have its file (gds_dir, gds_pattern).")

# ===========================================================================
# Step 1-2: load bundle + fit the null model
# ===========================================================================
if (!file.exists(pheno_path))
  stop(sprintf("Pheno bundle not found: %s\nRun 00-data-prep/assemble-pheno-covar.R first.", pheno_path))
pheno_data <- readRDS(pheno_path)
log_msg(sprintf("Bundle: %d samples (excluded %d / %d); covariates (%d): %s",
                length(pheno_data$Y), null_or(pheno_data$n_excluded, NA),
                null_or(pheno_data$n_total, NA), ncol(pheno_data$X),
                paste(pheno_data$covar_names, collapse = ", ")))

null_model <- fit_null_model(
  X         = pheno_data$X,
  Y         = pheno_data$Y,
  trait     = pheno_data$trait,
  sample_id = pheno_data$sample_id)
log_msg("Null model fitted.")

# Scan one chromosome.
scan_one <- function(k) {
  gds_path <- gds_path_for(k)
  if (!file.exists(gds_path))
    stop(sprintf("GDS file missing for chromosome %d: %s", k, gds_path))
  # The file must hold chromosome k: a mislabeled input stops here, before the scan.
  g <- SeqArray::seqOpen(gds_path, readonly = TRUE)
  chr_in_file <- unique(sub("^[Cc][Hh][Rr]", "", as.character(SeqArray::seqGetData(g, "chromosome"))))
  SeqArray::seqClose(g)
  if (!identical(chr_in_file, as.character(k)))
    stop(sprintf("%s holds chromosome(s) %s, not chromosome %d; check gds_files against chromosomes.",
                 gds_path, paste(chr_in_file, collapse = ", "), k))
  log_msg(sprintf("  Scanning chr%d...", k))
  t_chr <- proc.time()
  res <- marginal_scan(
    gds_file           = gds_path,
    null_model         = null_model,
    use_SPA            = use_SPA,
    chunk_size         = chunk_size,
    mac_cutoff         = mac_cutoff,
    missing_imputation = missing_imputation,
    output_csv         = file.path(results_dir, sprintf("marginal_chr%d.csv", k)),
    n_cores            = n_cores,
    verbose            = 1)
  n_var <- if (!is.null(res)) nrow(res) else 0L
  log_msg(sprintf("    chr%d: %d variants, %.1f s", k, n_var, (proc.time() - t_chr)[3]))
  invisible(res)
}

# ===========================================================================
# --chr N: scan only chr N (array-task unit; no combine / no plots)
# ===========================================================================
if (!is.null(chr)) {
  scan_one(chr)
  log_msg(sprintf("=== Marginal Scan (chr%d) Complete in %.1f s ===",
                  chr, (proc.time() - t_total_start)[3]))
  quit(save = "no", status = 0L)
}

# ===========================================================================
# Default: scan all `chroms`, then combine + plot.
# Bind the in-memory marginal_scan() returns, in chroms order.
# ===========================================================================
log_msg(sprintf("Scanning chromosomes: %s", paste(chroms, collapse = ", ")))
results_list <- vector("list", length(chroms))
for (j in seq_along(chroms)) results_list[[j]] <- scan_one(chroms[j])
combine_and_plot(in_memory = results_list)
log_msg(sprintf("=== Marginal Scan Complete in %.1f s (%.1f min) ===",
                (proc.time() - t_total_start)[3], (proc.time() - t_total_start)[3] / 60))
