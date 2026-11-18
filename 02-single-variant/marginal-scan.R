#!/usr/bin/env Rscript
# ============================================================================
# 02-single-variant/marginal-scan.R - per-variant association scan
#
# Cohort-agnostic per-variant (single-variant) association scan. Reads the
# assemble_pheno_covar bundle (pheno_path), fits the null model once
# (fit_null_model), then scans each chromosome's (a)GDS with GLOWr::marginal_scan
# (per-chr CSV) and, with --combine, rbinds them into single_variant_all.csv
# and draws the Manhattan/QQ plots.
#
# Modes (one of the two is given; each call is one unit of work):
#   --chr N      scan chromosome N -> results/single_variant_chr<N>.csv (the part
#                `scan`; no combine, no plots).
#   --combine    no scan: rbind the results/single_variant_chr*.csv of `chroms`
#                (numeric chr order) -> results/single_variant_all.csv + plots (the part
#                `scan_combine`).
#
# Output (under <output_dir> = <output_root>/<dataset>/<analysis>/single_variant/):
#   results/single_variant_chr<N>.csv   per-chromosome scan
#   results/single_variant_all.csv      combined (--combine)
#   plots/single_variant_plots.pdf      Manhattan + QQ (--combine)
#   logs/scan_log.txt             timestamped run log
#
# Usage (from the GLOWanalyses directory):
#   ./glow run single_variant --config <study>.yml --part scan     # the usual call: every unit of the part
#   Rscript 02-single-variant/marginal-scan.R --config <study>.yml --chr 22     # this script, for one unit
#   Rscript 02-single-variant/marginal-scan.R --config <study>.yml --combine     # this script, for one unit
# With --combine the script is the part scan_combine, which joins the per-chromosome tables.

suppressMessages(library(GLOWr))
library(SeqArray)
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <path> [--chr N | --combine] ----
source("lib/glow-stage.R")
st <- glow_stage_begin("single_variant", "scan", options = c(chr = "value", combine = "flag"),
                       part_if = c(combine = "scan_combine"))
config_path  <- st$config_path
chr          <- if (is.null(st$opts$chr)) NULL else as.integer(st$opts$chr)
combine_only <- isTRUE(st$opts$combine)
if (!is.null(chr) && combine_only) stop("--chr and --combine are mutually exclusive.")
if (!is.null(chr)) stopifnot(chr >= 1L, chr <= 22L)
# The stem of the result files, "single_variant" (release 0.3.6 wrote "marginal").
sv_stem      <- st$sv_stem
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

# ---- Run identity + output dirs (derived from the config path) ----
output_dir  <- st$output_dir
run_name    <- st$run_name
results_dir <- file.path(output_dir, "results")
plots_dir   <- file.path(output_dir, "plots")
logs_dir    <- file.path(output_dir, "logs")
for (d in c(results_dir, plots_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Required cohort fields + scan knobs (settings of the configuration) ----
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
# Combine helper: bind the per-chr tables (numeric chr order) -> single_variant_all.csv + plots.
# data.table's fread, rbindlist and fwrite (data.table is a GLOWr dependency), because a
# whole genome is about 125 million rows: read.csv, do.call(rbind) and write.csv took nine
# minutes for 11 million rows and hold several copies at once. The per-chromosome tables
# are released before the plots are drawn.
#
# The combined table is built from the per-chromosome FILES on every route. The default
# route writes them (marginal_scan(output_csv)) and reads them back here, exactly as the
# Terra route's Calibrate task and a cluster job array do with --combine, so the two routes
# write byte-identical tables. Binding the in-memory results instead gave values that differ
# from the round-tripped ones in the 15th significant digit (a few percent of the cells, a
# relative difference below 1e-14), which the route comparison of the image build rejects.
# Every chromosome of `chroms` must have its table: a missing one stops the combine, so an
# incomplete scatter cannot produce a genome with a chromosome silently absent.
# ===========================================================================
combine_and_plot <- function() {
  log_msg("Combining per-chromosome results")
  chr_files <- file.path(results_dir, sprintf("%s_chr%d.csv", sv_stem, chroms))
  absent    <- chr_files[!file.exists(chr_files)]
  if (length(absent) > 0L)
    stop("Per-chromosome table(s) missing for requested chromosome(s): ",
         paste(basename(absent), collapse = ", "),
         ". Every chromosome of `chroms` must be scanned before the combine.")
  results_list <- lapply(chr_files, function(f)
    data.table::fread(f, colClasses = list(character = "chr")))
  all_results <- data.table::rbindlist(results_list, use.names = TRUE)
  rm(results_list)                      # the per-chromosome copies go before the plots
  data.table::setDF(all_results)
  invisible(gc())

  csv_all <- file.path(results_dir, paste0(sv_stem, "_all.csv"))
  write_csv_exact(all_results, csv_all)
  log_msg(sprintf("  Combined: %d total variants -> %s", nrow(all_results), csv_all))

  # ---- Summary statistics ----
  log_msg(sprintf("  Variants with pvalue < 5e-8: %d",
                  sum(all_results$pvalue < 5e-8, na.rm = TRUE)))
  log_msg(sprintf("  Variants with pvalue < 1e-5: %d",
                  sum(all_results$pvalue < 1e-5, na.rm = TRUE)))

  # ---- Plots ----
  log_msg("Generating plots")
  pdf_path <- file.path(plots_dir, paste0(sv_stem, "_plots.pdf"))
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
  # the completion record of the part scan_combine
  glow_stage_end(st, outputs = c(file.path(results_dir, paste0(sv_stem, "_all.csv")),
                                 file.path(plots_dir, paste0(sv_stem, "_plots.pdf"))))
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
    output_csv         = file.path(results_dir, sprintf("%s_chr%d.csv", sv_stem, k)),
    n_cores            = n_cores,
    verbose            = 1)
  n_var <- if (!is.null(res)) nrow(res) else 0L
  log_msg(sprintf("    chr%d: %d variants, %.1f s", k, n_var, (proc.time() - t_chr)[3]))
  invisible(res)
}

# ===========================================================================
# --chr N: scan only chr N (one unit of work; no combine, no plots)
# ===========================================================================
if (!is.null(chr)) {
  scan_one(chr)
  log_msg(sprintf("=== Marginal Scan (chr%d) Complete in %.1f s ===",
                  chr, (proc.time() - t_total_start)[3]))
  # the completion record of this chromosome's unit of the part scan
  glow_stage_end(st, outputs = file.path(results_dir, sprintf("%s_chr%d.csv", sv_stem, chr)))
  quit(save = "no", status = 0L)
}

# A call without --chr and without --combine does not get here: the stage helper stops it,
# because the scan runs one chromosome per call (release 0.3.6 looped over all of them here).
stop("marginal-scan.R scans one chromosome per call (--chr N) or combines the tables (--combine).")
