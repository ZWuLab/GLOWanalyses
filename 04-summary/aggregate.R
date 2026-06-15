#!/usr/bin/env Rscript
# ============================================================================
# 04-summary/aggregate.R - Stage 3: region_type-aware aggregation
#
# Cohort-agnostic. Reads the per-unit Stage-2 flat tables and dispatches on the
# snapshot's region_type to the matching packaged aggregator, writing the
# aggregated tables under <output_dir>/aggregated/:
#   - gene   : read_all_chr_results() -> glow_results_all.csv + a top-K ranking
#              (glow_top_hits.csv with rank + passes_genome_wide).
#   - window : aggregate_scan_results() (genome rbind + alpha/n_windows threshold
#              + interval-union locus merge) -> scan_results_all.csv,
#              scan_sig_windows.csv, scan_loci.csv. Guarded by a completeness gate.
#   - coding : aggregate_coding_results() (3-source rbind + *_glowG rename +
#              (gene,category) merge + per-category significance) -> the three
#              source tables + coding_merged.csv + the two significance tables.
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 04-summary/aggregate.R \
#       --config runs/example/config.R \
#       [--top-k 100] [--allow-incomplete]

suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; top_k <- NA_integer_; allow_incomplete <- FALSE
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config"       && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--top-k"       && i < length(args)) { top_k <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--allow-incomplete") { allow_incomplete <- TRUE; i <- i + 1L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))

# ---- Setup (output_dir DERIVED; snapshot supplies config VALUES) ----
output_dir  <- file.path(dirname(config_path), "outputs")
run_name    <- basename(dirname(config_path))
cfg         <- load_config_snapshot(output_dir)
region_type <- match.arg(null_or(cfg$region_type, "gene"), c("gene", "window", "coding"))
out_dir     <- file.path(output_dir, "aggregated")
ts <- function() format(Sys.time(), "%H:%M:%S")
cat(sprintf("[%s] Aggregating run '%s' [region_type = %s]\n", ts(), run_name, region_type))

# ===========================================================================
if (region_type == "gene") {
  # --- Genome-wide (or subset-chr) rbind + top-K ranking ---
  # read_all_chr_results() takes the run's chromosomes (default 1:22); pass
  # cfg$chroms so a subset run aggregates exactly its chromosomes.
  results_dir <- file.path(output_dir, "results")
  all <- read_all_chr_results(results_dir, chroms = sort(null_or(cfg$chroms, 1:22)))
  n_results <- nrow(all)
  if (n_results == 0L) stop("No per-chr gene results in ", file.path(output_dir, "results"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(all, file.path(out_dir, "glow_results_all.csv"), row.names = FALSE)

  primary   <- null_or(cfg$primary_test, "GLOW_Omni")
  top_k     <- if (is.na(top_k)) null_or(cfg$top_k, 100L) else top_k
  threshold <- null_or(cfg$alpha, 0.05) / n_results
  top <- all[order(all[[primary]]), , drop = FALSE]
  top$rank <- seq_len(nrow(top))
  top$passes_genome_wide <- top[[primary]] < threshold
  top <- head(top, top_k)
  write.csv(top, file.path(out_dir, "glow_top_hits.csv"), row.names = FALSE)
  cat(sprintf("[%s] %d genes; Bonferroni (%.3g/%d) = %.3g; %d of top-%d pass.\n",
              ts(), n_results, null_or(cfg$alpha, 0.05), n_results, threshold,
              sum(top$passes_genome_wide), nrow(top)))

# ===========================================================================
} else if (region_type == "window") {
  output_format <- null_or(cfg$output_format, "csv")
  ext <- if (identical(output_format, "fst")) "fst" else "csv"
  # --- Completeness gate (a missing chunk deflates the alpha/n_windows denom) ---
  comp <- check_scan_completeness(output_dir, output_format)
  cat(sprintf("[%s] completeness: %d/%d chunks present (%d missing); status %s\n",
              ts(), comp$n_present, comp$n_expected, length(comp$missing_ids),
              if (comp$complete) "COMPLETE" else "INCOMPLETE"))
  if (!comp$complete && !allow_incomplete)
    stop(sprintf("INCOMPLETE scan: %d/%d chunks. Resubmit missing chunks or pass ",
                 comp$n_present, comp$n_expected), "--allow-incomplete.")

  chunk_files <- list.files(file.path(output_dir, "results"),
                            pattern = sprintf("^scan_chunk_\\d+\\.%s$", ext),
                            full.names = TRUE)
  if (length(chunk_files) == 0L)
    stop("No per-chunk tables (scan_chunk_*.", ext, ") in ", file.path(output_dir, "results"))
  agg <- aggregate_scan_results(
    chunk_files,
    primary_test  = null_or(cfg$primary_test, "GLOW_Omni"),
    alpha         = null_or(cfg$alpha, 0.05),
    merge_gap     = null_or(cfg$merge_gap, 0L),
    cmac_cutoff   = cfg$cmac_cutoff,
    output_format = output_format)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(agg$scan_results_all, file.path(out_dir, "scan_results_all.csv"), row.names = FALSE)
  write.csv(agg$scan_sig_windows, file.path(out_dir, "scan_sig_windows.csv"), row.names = FALSE)
  write.csv(agg$scan_loci,        file.path(out_dir, "scan_loci.csv"),        row.names = FALSE)
  s <- agg$summary
  cat(sprintf("[%s] n_windows = %d; threshold = %.3g; %d sig window(s) -> %d locus/loci (%s).\n",
              ts(), s$n_windows, s$threshold, s$n_sig_windows, s$n_loci, s$primary_test))

# ===========================================================================
} else if (region_type == "coding") {
  staar_modes  <- null_or(cfg$staar_modes, c("spa", "nospa"))
  primary_test <- null_or(cfg$primary_test_glow, "GLOW_Omni")
  # aggregate_coding_results owns the 3-source rbind, the *_glowG rename, the
  # (gene,category) merge, the per-category significance, and the completeness
  # reconciliation. Run once; write every table it returns.
  agg <- aggregate_coding_results(
    output_dir        = output_dir,
    chroms            = cfg$chroms,
    gene_num_in_array = cfg$gene_num_in_array,
    staar_modes       = staar_modes,
    alpha             = null_or(cfg$alpha, 0.05),
    gw_alpha          = null_or(cfg$staar_genomewide_alpha, 2.5e-6),
    primary_test      = primary_test,
    allow_partial     = allow_incomplete)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Source tables.
  write.csv(agg$glow_coding_all, file.path(out_dir, "glow_coding_all.csv"), row.names = FALSE)
  for (mode in staar_modes) {
    flat <- agg[[paste0("staar_coding_", mode, "_all")]]
    if (!is.null(flat))
      write.csv(flat, file.path(out_dir, sprintf("staar_coding_%s_all.csv", mode)), row.names = FALSE)
  }
  # Merge + significance tables.
  write.csv(agg$coding_merged, file.path(out_dir, "coding_merged.csv"), row.names = FALSE)
  write.csv(agg$coding_sig_by_category, file.path(out_dir, "coding_sig_by_category.csv"), row.names = FALSE)
  write.csv(agg$coding_sig_by_category_cells, file.path(out_dir, "coding_sig_by_category_cells.csv"), row.names = FALSE)
  cat(sprintf("[%s] GLOW: %d (gene,category) cells; merged %d cells; %d (category,method) sig rows; %d sig cells.\n",
              ts(), nrow(agg$glow_coding_all), nrow(agg$coding_merged),
              nrow(agg$coding_sig_by_category), nrow(agg$coding_sig_by_category_cells)))
}

cat(sprintf("[%s] Aggregation complete. Tables in: %s\n", ts(), out_dir))
