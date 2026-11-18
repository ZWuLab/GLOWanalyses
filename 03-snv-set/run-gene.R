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
# Output (under <output_dir>/results/): glow_chr<N>.csv. With staar_enabled = TRUE in
# the config the table also carries the 8 STAAR comparator columns, and
# write_staar_detail = TRUE adds the staar_detail/ sidecar.
#
# Usage (from the GLOWanalyses directory):
#   ./glow run snv_set --config <study>.yml --part scan     # the usual call: every unit of the part
#   Rscript 03-snv-set/run-gene.R --config <study>.yml --region-type gene --chr 22     # this script, for one unit
# A value that release 0.3.6 took as an option here is a setting of the configuration: --max-genes is snv_set.max_genes, --n-cores is shared.n_cores.

suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
library(SeqArray)
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <path> --chr <N> [--max-genes K] [--n-cores N] ----
source("lib/glow-stage.R")
st <- glow_stage_begin("snv_set", "scan", load = FALSE,
                       options = c(chr = "value", "max-genes" = "value", "n-cores" = "value"))
config_path <- st$config_path
chr       <- if (is.null(st$opts$chr)) NULL else as.integer(st$opts$chr)
max_genes <- as.integer(glow_stage_option(st, "max-genes", c("snv_set", "max_genes"), NA_integer_))
# The workers of this unit: the current `shared.n_cores` of the configuration, so that a
# changed number reaches new units without a new preparation.
n_cores   <- glow_stage_option(st, "n-cores", c("shared", "n_cores"), NULL)
if (!is.null(n_cores)) n_cores <- as.integer(n_cores)
stopifnot(!is.null(chr), chr >= 1L, chr <= 22L)

# ---- Run identity + Stage-1 artifacts ----
output_dir <- st$output_dir
run_name    <- st$run_name
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
if (!is.null(n_cores))
  log_msg(sprintf("n_cores = %d from the current setting or --n-cores (the Stage-1 snapshot says %s)",
                  n_cores, format(null_or(config$n_cores, 1L))))

# ---- Models + per-chr inputs ----
pi_models <- load_PI_models(config$pi_model_dir)
b_model   <- if (!is.null(config$b_model_path)) load_B_model(config$b_model_path) else NULL

# ---- Optional embedded STAAR comparator (non-SPA STAAR::STAAR on GLOW's exact G) ----
# Stage 1 fitted 01-staar_null_model.rds on the same phenotype and covariates as
# GLOW's null model when the config set staar_enabled <- TRUE. The context makes
# run_scan_unit() fill the 8 STAAR family columns for every gene (all NA otherwise).
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
# What an earlier run of this unit left goes first, so that the
# completion record lists the files of this run.
glow_scan_unit_files(st, results_dir, sprintf("glow_chr%d", chr), sprintf("chr%d", chr), clear = TRUE)
run_scan_unit(
  gds_path = gds_path, null_model = null_model, pi_models = pi_models,
  pi_features = config$pi_features, reference_medians = ref_medians,
  b_func = config$b_func, b_model = b_model,
  region_type = "gene", config = config,
  results_path = file.path(results_dir, sprintf("glow_chr%d.csv", chr)),
  output_dir = output_dir, output_policy = policy,
  chr = chr, unit_id = sprintf("chr%d", chr),
  staar_context = staar_context, use_spa = config$use_spa, z_scale = z_scale,
  debug_errors = isTRUE(config$debug_errors),
  max_regions = max_genes, n_cores = n_cores, progress_every = 50L, logger = log_msg, verbose = 0)

log_msg(sprintf("Done gene chr%d [run %s].", chr, run_name))

# The completion record of this unit.
glow_stage_end(st, outputs = glow_scan_unit_files(st, results_dir, sprintf("glow_chr%d", chr), sprintf("chr%d", chr)))
