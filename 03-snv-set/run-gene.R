#!/usr/bin/env Rscript
# ============================================================================
# 03-snv-set/run-gene.R - Stage 2 (gene), per chromosome
#
# Thin entry over the unified Stage-2 driver: parse args, load the Stage-1
# snapshot, compute this chromosome's annotation medians (gene carries no
# Stage-1 median cache), then call GLOWpipeline::run_scan_unit(region_type =
# "gene") - the SAME scan the window/coding runs use. Writes the per-chr
# gene-level flat table + policy-gated sidecars + the per-chr log.
#
# Output (under <output_dir>/results/): glow_chr<N>.csv
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 03-snv-set/run-gene.R \
#       --config runs/example/config.R \
#       --chr 22 [--max-genes 20]

suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
library(SeqArray)
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <path> --chr <N> [--max-genes K] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; chr <- NULL; max_genes <- NA_integer_
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config"      && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--chr"        && i < length(args)) { chr <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--max-genes"  && i < length(args)) { max_genes <- as.integer(args[i + 1L]); i <- i + 2L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path), !is.null(chr), chr >= 1L, chr <= 22L)

# ---- Run identity + Stage-1 artifacts ----
output_dir  <- file.path(dirname(config_path), "outputs")
run_name    <- basename(dirname(config_path))
shared_dir  <- file.path(output_dir, "shared")
results_dir <- file.path(output_dir, "results")
logs_dir    <- file.path(output_dir, "logs")
for (d in c(results_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

need <- file.path(shared_dir, c("01-config_used.rds", "01-null_model.rds"))
if (any(!file.exists(need)))
  stop("Missing Stage-1 artifacts. Run 03-snv-set/prepare.R first. Missing:\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "))
config     <- readRDS(file.path(shared_dir, "01-config_used.rds"))
null_model <- readRDS(file.path(shared_dir, "01-null_model.rds"))

# ---- Logging ----
log_con <- file(file.path(logs_dir, sprintf("log_chr%d.txt", chr)), open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE, after = FALSE)
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = ""); writeLines(line, log_con); flush(log_con)
}
log_msg(sprintf("Starting gene chr%d [run %s]; n_samples = %d, trait = %s",
                chr, run_name, null_model$n, null_model$trait))

# ---- Models + per-chr inputs ----
pi_models <- load_PI_models(config$pi_model_dir)
b_model   <- if (!is.null(config$b_model_path)) load_B_model(config$b_model_path) else NULL
gds_path  <- file.path(config$gds_dir, gsub("\\{chr\\}", chr, config$gds_pattern))
stopifnot(file.exists(gds_path))

# Gene carries no Stage-1 median cache: compute this chromosome's medians here.
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

# ---- Run the unified Stage-2 driver (region_type = "gene") ----
run_scan_unit(
  gds_path = gds_path, null_model = null_model, pi_models = pi_models,
  pi_features = config$pi_features, reference_medians = ref_medians,
  b_func = config$b_func, b_model = b_model,
  region_type = "gene", config = config,
  results_path = file.path(results_dir, sprintf("glow_chr%d.csv", chr)),
  output_dir = output_dir, output_policy = policy,
  chr = chr, unit_id = sprintf("chr%d", chr),
  use_spa = config$use_spa, z_scale = z_scale,
  max_regions = max_genes, progress_every = 50L, logger = log_msg, verbose = 0)

log_msg(sprintf("Done gene chr%d [run %s].", chr, run_name))
