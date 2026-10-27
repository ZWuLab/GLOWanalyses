#!/usr/bin/env Rscript
# ============================================================================
# 03-snv-set/run-noncoding.R - Stage 2 (noncoding), per chromosome
#
# Thin entry over the unified Stage-2 driver for the gene-centric non-coding
# masks: parse args, load the Stage-1 snapshot (and the STAAR null model when
# staar_enabled), compute this chromosome's annotation medians, then call
# GLOWpipeline::run_scan_unit(region_type = "noncoding"), the SAME scan the
# gene/window/coding runs use. run_scan_unit() builds one regions table per
# mask of config$noncoding_categories from the chromosome's annotations
# (GLOWr::define_regions_noncoding; the five fields genecode_comprehensive_category,
# genecode_comprehensive_info, genehancer, cage_tc and rdhs), loops the masks
# as the coding run loops its categories, tags each row with `category`, and
# writes the per-(gene, mask) flat table in gene-first order.
#
# Output (under <output_dir>/results-glow/): glow_noncoding_chr<N>.csv, and
# glow_noncoding_chr<N>_errors.csv when a cell errored. With staar_enabled =
# TRUE in the config the table also carries the 8 STAAR comparator columns on
# GLOW's own set of every cell.
#
# Usage (from the GLOWanalyses directory):
#   conda activate r_env
#   Rscript 03-snv-set/run-noncoding.R \
#       --config runs/<name>/config.R \
#       --chr 22 [--max-genes 20] [--n-cores N]   # --max-genes caps every mask's table;
#                                                  # --n-cores overrides the snapshot's n_cores

suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
library(SeqArray)
null_or <- function(a, b) if (is.null(a)) b else a
t_start <- Sys.time()

# ---- Args: --config <path> --chr <N> [--max-genes K] [--n-cores N] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; chr <- NULL; max_genes <- NA_integer_; n_cores <- NULL
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config"      && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--chr"        && i < length(args)) { chr <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--max-genes"  && i < length(args)) { max_genes <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--n-cores"    && i < length(args)) { n_cores <- as.integer(args[i + 1L]); i <- i + 2L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path), !is.null(chr), chr >= 1L, chr <= 22L)

# ---- Run identity + Stage-1 artifacts ----
output_dir       <- file.path(dirname(config_path), "outputs")
run_name         <- basename(dirname(config_path))
shared_dir       <- file.path(output_dir, "shared")
results_glow_dir <- file.path(output_dir, "results-glow")
logs_dir         <- file.path(output_dir, "logs")
for (d in c(results_glow_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

need <- file.path(shared_dir, c("01-config_used.rds", "01-null_model.rds"))
if (any(!file.exists(need)))
  stop("Missing Stage-1 artifacts. Run 03-snv-set/prepare.R first. Missing:\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "))
config     <- readRDS(file.path(shared_dir, "01-config_used.rds"))
null_model <- readRDS(file.path(shared_dir, "01-null_model.rds"))
if (!identical(config$region_type, "noncoding"))
  stop("Stage 1 was prepared for region_type '", config$region_type, "', not 'noncoding'.")
# The per-chromosome table is named .csv here and read back by read.csv in
# 04-summary/aggregate.R, so the fst format cannot be aggregated on this route.
if (!identical(null_or(config$output_format, "csv"), "csv"))
  stop("region_type 'noncoding' writes and aggregates CSV tables only; set output_format <- ",
       "\"csv\" in the run config (it asks for '", config$output_format, "').")

# ---- Logging ----
log_con <- file(file.path(logs_dir, sprintf("log_noncoding_chr%d.txt", chr)), open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE, after = FALSE)
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = ""); writeLines(line, log_con); flush(log_con)
}
masks <- null_or(config$noncoding_categories, noncoding_masks())
log_msg(sprintf("Starting noncoding chr%d [run %s]: %d masks (%s); n_samples = %d, trait = %s",
                chr, run_name, length(masks), paste(masks, collapse = ","),
                null_model$n, null_model$trait))
if (!is.null(n_cores))
  log_msg(sprintf("n_cores = %d from --n-cores (the Stage-1 snapshot says %s)",
                  n_cores, format(null_or(config$n_cores, 1L))))

# ---- Models + per-chr inputs ----
pi_models <- load_PI_models(config$pi_model_dir)
b_model   <- if (!is.null(config$b_model_path)) load_B_model(config$b_model_path) else NULL

# ---- Optional embedded STAAR comparator (non-SPA STAAR::STAAR on GLOW's exact G) ----
# Stage 1 fitted 01-staar_null_model.rds on the same phenotype and covariates as
# GLOW's null model when the config set staar_enabled <- TRUE. The context makes
# run_scan_unit() fill the 8 STAAR family columns for every (gene, mask) cell.
staar_context <- NULL
if (isTRUE(config$staar_enabled)) {
  staar_null_path <- file.path(shared_dir, "01-staar_null_model.rds")
  if (!file.exists(staar_null_path))
    stop("Stage-1 says staar_enabled but no STAAR null model at:\n  ", staar_null_path)
  suppressMessages(library(STAAR))
  staar_context <- list(
    null_model      = readRDS(staar_null_path),
    anno_cols       = null_or(config$staar_anno_features, config$pi_features),
    rare_maf_cutoff = config$filter_spec$rare_maf_cutoff)
  log_msg(sprintf("STAAR comparator ON (non-SPA; %d annotation features; same covariates as GLOW)",
                  length(staar_context$anno_cols)))
}
gds_path <- file.path(config$gds_dir, gsub("\\{chr\\}", chr, config$gds_pattern))
stopifnot(file.exists(gds_path))

# Noncoding carries no Stage-1 median cache: compute this chromosome's medians
# here, under the run's filter, as the gene run does.
filter <- do.call(variant_filter, config$filter_spec)
log_msg("Computing chromosome-wide annotation medians...")
ref_medians <- compute_annotation_medians(
  gds_path, config$pi_features, filter, sample_id = null_model$sample_id, verbose = 0)

# ---- Output policy + GC calibration ----
policy <- glow_output_policy(
  write_flat_table     = TRUE,
  write_evidence       = isTRUE(config$write_evidence),
  write_staar_detail   = isTRUE(config$write_staar_detail),
  write_result_objects = isTRUE(config$write_result_objects))
inflation_factor <- null_or(config$inflation_factor, 1.0)
z_scale <- 1 / sqrt(inflation_factor)
if (inflation_factor != 1.0)
  log_msg(sprintf("GC calibration: inflation_factor = %.4f -> z_scale = %.4f", inflation_factor, z_scale))

# ---- Run the unified Stage-2 driver (region_type = "noncoding") ----
run_scan_unit(
  gds_path = gds_path, null_model = null_model, pi_models = pi_models,
  pi_features = config$pi_features, reference_medians = ref_medians,
  b_func = config$b_func, b_model = b_model,
  region_type = "noncoding", config = config,
  results_path = file.path(results_glow_dir, sprintf("glow_noncoding_chr%d.csv", chr)),
  output_dir = output_dir, output_policy = policy,
  chr = chr, unit_id = sprintf("chr%d", chr),
  staar_context = staar_context, use_spa = config$use_spa, z_scale = z_scale,
  max_regions = max_genes, n_cores = n_cores, progress_every = 200L, logger = log_msg, verbose = 0)

log_msg(sprintf("Done noncoding chr%d [run %s]. masks=%d, wall=%.1fs",
                chr, run_name, length(masks),
                as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
