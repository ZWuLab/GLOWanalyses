#!/usr/bin/env Rscript
# ============================================================================
# 04-summary/plots.R - Stage 3: Manhattan / QQ + per-test lambda (+ drilldown)
#
# Cohort-agnostic. Reads the aggregated table from 04-summary/aggregate.R and,
# dispatching on the snapshot's region_type, writes the headline visuals under
# <output_dir>/aggregated/:
#   - gene/window : a Manhattan of the primary test (genome-wide Bonferroni line,
#       top hits highlighted) + the per-test lambda + QQ pass via
#       GLOWpipeline::summarize_qq_and_lambdas() (families via
#       GLOWr::group_glow_tests_by_family()). For gene, when evidence packets
#       were written, a per-top-hit SNV-level drilldown via GLOWr::glow_drilldown_df().
#   - coding : per-category lambda + a 4-method QQ overlay on the >=2-variant
#       common support, via summarize_qq_and_lambdas() (methods as one "family").
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 04-summary/plots.R \
#       --config runs/example/config.R [--drilldown]

suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
null_or <- function(a, b) if (is.null(a)) b else a
ts <- function() format(Sys.time(), "%H:%M:%S")

# ---- Args ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; do_drilldown <- FALSE
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--drilldown") { do_drilldown <- TRUE; i <- i + 1L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))

# ---- Setup ----
output_dir  <- file.path(dirname(config_path), "outputs")
cfg         <- load_config_snapshot(output_dir)
region_type <- match.arg(null_or(cfg$region_type, "gene"), c("gene", "window", "coding"))
out_dir     <- file.path(output_dir, "aggregated")
stopifnot(dir.exists(out_dir))
cat(sprintf("[%s] Plots for run '%s' [region_type = %s]\n", ts(), cfg$run_name, region_type))

# ===========================================================================
if (region_type %in% c("gene", "window")) {
  primary <- null_or(cfg$primary_test, "GLOW_Omni")
  alpha   <- null_or(cfg$alpha, 0.05)
  results_csv <- if (region_type == "gene")
    file.path(out_dir, "glow_results_all.csv") else file.path(out_dir, "scan_results_all.csv")
  if (!file.exists(results_csv)) stop("Run 04-summary/aggregate.R first (missing ", results_csv, ").")

  df <- read.csv(results_csv, stringsAsFactors = FALSE, check.names = FALSE,
                 colClasses = c(chr = "character"))
  df$chr_int <- suppressWarnings(as.integer(df$chr))
  # x-axis: gene uses the region midpoint, window uses the always-on `midpoint`.
  df$pos <- if (region_type == "window") as.integer(df$midpoint)
            else as.integer((df$start + df$end) / 2L)
  df$variant_id <- df$label

  n_results  <- nrow(df)
  threshold  <- alpha / n_results
  top_labels <- select_highlight_regions(df, primary, threshold, min_n = 10L)
  n_pass     <- sum(df[[primary]] <= threshold, na.rm = TRUE)
  cat(sprintf("[%s] n = %d; Bonferroni = %.3g; %d pass; highlighting %d.\n",
              ts(), n_results, threshold, n_pass, length(top_labels)))

  # --- Manhattan ---
  pdf(file.path(out_dir, "manhattan_omni.pdf"), width = 12, height = 5)
  plot_manhattan(df, pvalue_col = primary, chr_col = "chr_int", pos_col = "pos",
                 suggestive_line = NULL, genome_wide_line = threshold,
                 highlight = top_labels,
                 title = sprintf("GLOW %s (%s trait)", region_type, null_or(cfg$trait, "?")))
  dev.off()

  # --- Per-test lambda + QQ (families auto-derived from the GLOW_* columns) ---
  glow_tests <- grep("^GLOW_", names(df), value = TRUE)
  families   <- group_glow_tests_by_family(glow_tests)
  res_all <- summarize_qq_and_lambdas(
    df, label = "all", output_dir = out_dir, omni_test = primary, families = families,
    subtitle = sprintf("all units, n = %d", n_results))
  cat(sprintf("[%s] lambda_%s = %.3f; per-test lambdas in .family_lambdas_all.rds.\n",
              ts(), primary, res_all$lambda_omni))

  # Window: also the multi-variant common support (apples-to-apples vs STAAR).
  if (region_type == "window" && "is_single" %in% names(df)) {
    df_mv <- df[!df$is_single, , drop = FALSE]
    summarize_qq_and_lambdas(
      df_mv, label = "multivariant", output_dir = out_dir, omni_test = primary,
      families = families,
      subtitle = sprintf("multi-variant units only, n = %d of %d", nrow(df_mv), n_results))
    cat(sprintf("[%s] multi-variant QQ/lambda also written (n = %d).\n", ts(), nrow(df_mv)))
  }

  # --- Per-top-hit SNV-level drilldown (gene; needs evidence packets) ---
  if (region_type == "gene" && do_drilldown) {
    ev_dir <- file.path(output_dir, "evidence")
    top_csv <- file.path(out_dir, "glow_top_hits.csv")
    if (!dir.exists(ev_dir)) {
      cat(sprintf("[%s] --drilldown skipped: no evidence/ (set write_evidence = TRUE).\n", ts()))
    } else if (!file.exists(top_csv)) {
      cat(sprintf("[%s] --drilldown skipped: no glow_top_hits.csv.\n", ts()))
    } else {
      top <- read.csv(top_csv, stringsAsFactors = FALSE)
      needed_chrs <- sort(unique(as.integer(top$chr)))
      evidence <- load_evidence_for_chrs(ev_dir, needed_chrs)
      dd_dir <- file.path(out_dir, "top_hits"); dir.create(dd_dir, showWarnings = FALSE)
      n_done <- 0L
      for (g in seq_len(nrow(top))) {
        e <- evidence[[top$label[g]]]
        if (is.null(e)) next
        write.csv(glow_drilldown_df(e),
                  file.path(dd_dir, sprintf("%s-variants.csv", top$label[g])), row.names = FALSE)
        n_done <- n_done + 1L
      }
      cat(sprintf("[%s] drilldown: %d/%d top-hit gene CSVs -> %s\n", ts(), n_done, nrow(top), dd_dir))
    }
  }

# ===========================================================================
} else if (region_type == "coding") {
  merged_path <- file.path(out_dir, "coding_merged.csv")
  if (!file.exists(merged_path)) stop("Run 04-summary/aggregate.R first (missing ", merged_path, ").")
  merged <- read.csv(merged_path, stringsAsFactors = FALSE, colClasses = c(chr = "character"))

  methods <- intersect(c("GLOW_Omni", "STAAR_O_glowG", "STAAR_O_native", "STAAR_B_native"),
                       names(merged))
  primary <- null_or(cfg$primary_test_glow, "GLOW_Omni")
  # Per (category, method) lambda on the >=2-variant common support (where all
  # methods are defined); the synonymous category is the calibration readout.
  common <- if ("is_single" %in% names(merged))
    merged[!is.na(merged$is_single) & !merged$is_single, , drop = FALSE] else merged
  cats <- sort(unique(merged$category))
  lam_rows <- list()
  for (cat0 in cats) {
    sc <- common[common$category == cat0, , drop = FALSE]
    for (meth in methods) {
      p <- sc[[meth]]
      lam_rows[[length(lam_rows) + 1L]] <- data.frame(
        category = cat0, method = meth,
        n = sum(is.finite(p) & p > 0 & p <= 1),
        lambda_gc = compute_lambda_gc(p), stringsAsFactors = FALSE)
    }
    # Per-category 4-method QQ overlay via summarize_qq_and_lambdas (methods as one family).
    if (length(methods) > 0L && nrow(sc) > 0L) {
      summarize_qq_and_lambdas(
        sc, label = cat0, output_dir = out_dir, omni_test = primary,
        families = setNames(list(methods), cat0),
        family_title = function(fam) sprintf("QQ %s", fam),
        pdf_basename_family = "qq", write_rds = FALSE,
        verbose = FALSE)
    }
  }
  lam <- do.call(rbind, lam_rows)
  write.csv(lam, file.path(out_dir, "lambda_by_category.csv"), row.names = FALSE)
  cat(sprintf("[%s] lambda_by_category.csv (%d rows) + per-category qq_<cat>.pdf written.\n",
              ts(), nrow(lam)))
}

cat(sprintf("[%s] Plots complete. Artifacts in: %s\n", ts(), out_dir))
