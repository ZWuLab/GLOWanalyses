#!/usr/bin/env Rscript
# ============================================================================
# 03-snv-set/run-coding.R - Stage 2 (coding), per chromosome
#
# Thin entry over the unified Stage-2 driver: parse args, load the Stage-1
# snapshot + the non-SPA STAAR null model (for the baked-in STAAR_O_glowG), then
# call GLOWpipeline::run_scan_unit(region_type = "coding") - the SAME scan the
# gene/window runs use. run_scan_unit() loops the configured coding categories
# OVER the chromosome's gene table (category-outer: one glow_scan_chunk call per
# category; the category -> filter map lives in build_scan_regions), tags each
# row with `category`, rbinds, and writes the per-(gene, category) flat table in
# gene-first order.
#
# Output (under <output_dir>/results-glow/): glow_coding_chr<N>.csv
#
# A `noncoding` slot (gene-centric non-coding RNA / regulatory sets) is FUTURE
# work: add one more entry script (run-noncoding.R) + a build_scan_regions
# dispatcher case, NOT new machinery.
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 03-snv-set/run-coding.R \
#       --config runs/example/config.R \
#       --chr 22 [--max-genes 20]

suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
library(SeqArray)
null_or <- function(a, b) if (is.null(a)) b else a
t_start <- Sys.time()

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
output_dir       <- file.path(dirname(config_path), "outputs")
run_name         <- basename(dirname(config_path))
shared_dir       <- file.path(output_dir, "shared")
results_glow_dir <- file.path(output_dir, "results-glow")
logs_dir         <- file.path(output_dir, "logs")
for (d in c(results_glow_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

need <- file.path(shared_dir, c("01-config_used.rds", "01-null_model.rds",
                                "01-annotation_medians.rds"))
if (any(!file.exists(need)))
  stop("Missing Stage-1 artifacts. Run 03-snv-set/prepare.R first. Missing:\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "))
config             <- readRDS(file.path(shared_dir, "01-config_used.rds"))
null_model         <- readRDS(file.path(shared_dir, "01-null_model.rds"))
ref_medians_by_chr <- readRDS(file.path(shared_dir, "01-annotation_medians.rds"))

# ---- Logging ----
log_con <- file(file.path(logs_dir, sprintf("log_coding_chr%d.txt", chr)), open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE, after = FALSE)
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = ""); writeLines(line, log_con); flush(log_con)
}

# ---- Median cache + GDS for this chromosome ----
ref_med <- ref_medians_by_chr[[as.character(chr)]]
if (is.null(ref_med))
  stop(sprintf("No cached medians for chr%d. Re-run prepare.R with this chromosome.", chr))
gds_path <- file.path(config$gds_dir, gsub("\\{chr\\}", chr, config$gds_pattern))
stopifnot(file.exists(gds_path))

category_specs  <- config$categories
category_labels <- if (is.list(category_specs)) names(category_specs) else category_specs
log_msg(sprintf("Starting coding chr%d [run %s]: %d categories (%s)",
                chr, run_name, length(category_labels), paste(category_labels, collapse = ",")))

# ---- Models + output policy ----
pi_models <- load_PI_models(config$pi_model_dir)
b_model   <- if (!is.null(config$b_model_path)) load_B_model(config$b_model_path) else NULL
policy <- glow_output_policy(
  write_flat_table     = TRUE,
  write_evidence       = isTRUE(config$write_evidence),
  write_staar_detail   = isTRUE(config$write_staar_detail),
  write_result_objects = isTRUE(config$write_result_objects))

# ---- Baked-in STAAR_O_glowG context (non-SPA STAAR::STAAR on GLOW's exact G) ----
staar_glowG_ctx <- NULL
if (isTRUE(config$staar_bakein_glow)) {
  staar_null_path <- file.path(shared_dir, "01-staar_null_model_nospa.rds")
  if (!file.exists(staar_null_path))
    stop("staar_bakein_glow = TRUE but no non-SPA STAAR null model at:\n  ", staar_null_path)
  suppressMessages(library(STAAR))
  staar_glowG_ctx <- list(
    null_model      = readRDS(staar_null_path),
    anno_cols       = null_or(config$staar_anno_features, config$pi_features),
    rare_maf_cutoff = config$filter_spec$rare_maf_cutoff)
  log_msg(sprintf("Baked-in STAAR_O_glowG ON (non-SPA; %d anno features)",
                  length(staar_glowG_ctx$anno_cols)))
}

# ---- GC calibration + SPA option (GLOW-only; STAAR untouched) ----
inflation_factor <- null_or(config$inflation_factor, 1.0)
z_scale <- 1 / sqrt(inflation_factor)
if (inflation_factor != 1.0)
  log_msg(sprintf("GC calibration: inflation_factor = %.4f -> z_scale = %.4f (GLOW-only)",
                  inflation_factor, z_scale))
use_spa <- config$use_spa
if (!is.null(use_spa)) log_msg(sprintf("SPA option: use_spa = %s (GLOW-only)", use_spa))

# ---- Run the unified Stage-2 driver (region_type = "coding") ----
run_scan_unit(
  gds_path = gds_path, null_model = null_model, pi_models = pi_models,
  pi_features = config$pi_features, reference_medians = ref_med,
  b_func = config$b_func, b_model = b_model,
  region_type = "coding", config = config,
  results_path = file.path(results_glow_dir, sprintf("glow_coding_chr%d.csv", chr)),
  output_dir = output_dir, output_policy = policy,
  chr = chr, unit_id = sprintf("chr%d", chr),
  staar_context = staar_glowG_ctx, use_spa = use_spa, z_scale = z_scale,
  max_regions = max_genes, progress_every = 200L, logger = log_msg, verbose = 0)

log_msg(sprintf("Done coding chr%d [run %s]. categories=%d, wall=%.1fs",
                chr, run_name, length(category_labels),
                as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
