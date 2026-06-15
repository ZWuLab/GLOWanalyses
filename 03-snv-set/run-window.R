#!/usr/bin/env Rscript
# ============================================================================
# 03-snv-set/run-window.R - Stage 2 (window), per genome chunk
#
# Thin entry over the unified Stage-2 driver: parse args, load the Stage-1
# snapshot + median cache + chunk table, resolve THIS chunk + its inputs, then
# call GLOWpipeline::run_scan_unit(region_type = "window") - the SAME scan the
# gene/coding runs use. run_scan_unit() builds the chunk's window table (the
# boundary rule lives in build_scan_regions), runs the scan, and writes the
# per-chunk flat table + policy-gated sidecars + the per-chunk log.
#
# Output (under <output_dir>/results/): scan_chunk_<id>.csv (or .fst)
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 03-snv-set/run-window.R \
#       --config runs/example/config.R \
#       --chunk 1 [--max-windows 5]

suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
library(SeqArray)
null_or <- function(a, b) if (is.null(a)) b else a
t_start <- Sys.time()

# ---- Args: --config <path> --chunk <id> [--max-windows K] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; chunk_id <- NULL; max_windows <- NA_integer_
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config"       && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--chunk"       && i < length(args)) { chunk_id <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--max-windows" && i < length(args)) { max_windows <- as.integer(args[i + 1L]); i <- i + 2L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path), !is.null(chunk_id), chunk_id >= 1L)

# ---- Run identity + Stage-1 artifacts ----
output_dir  <- file.path(dirname(config_path), "outputs")
run_name    <- basename(dirname(config_path))
shared_dir  <- file.path(output_dir, "shared")
results_dir <- file.path(output_dir, "results")
logs_dir    <- file.path(output_dir, "logs")
for (d in c(results_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

need <- file.path(shared_dir, c("01-config_used.rds", "01-null_model.rds",
                                "01-annotation_medians.rds", "01-genome_chunks.rds"))
if (any(!file.exists(need)))
  stop("Missing Stage-1 artifacts. Run 03-snv-set/prepare.R first. Missing:\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "))
config             <- readRDS(file.path(shared_dir, "01-config_used.rds"))
null_model         <- readRDS(file.path(shared_dir, "01-null_model.rds"))
ref_medians_by_chr <- readRDS(file.path(shared_dir, "01-annotation_medians.rds"))
genome_chunks      <- readRDS(file.path(shared_dir, "01-genome_chunks.rds"))

# ---- Logging ----
log_con <- file(file.path(logs_dir, sprintf("log_chunk_%04d.txt", chunk_id)), open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE, after = FALSE)
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = ""); writeLines(line, log_con); flush(log_con)
}

# ---- Resolve THIS chunk ----
chunk <- genome_chunks[genome_chunks$chunk_id == chunk_id, , drop = FALSE]
if (nrow(chunk) != 1L)
  stop(sprintf("chunk_id %d not found in the chunk table (%d chunks).",
               chunk_id, nrow(genome_chunks)))
chr <- as.character(chunk$chr)
log_msg(sprintf("Starting chunk %d [run %s]: chr%s [%d, %d], est_n_windows=%d",
                chunk_id, run_name, chr, chunk$chunk_start, chunk$chunk_end, chunk$est_n_windows))

ref_medians <- ref_medians_by_chr[[chr]]
if (is.null(ref_medians))
  stop(sprintf("No cached medians for chr%s. Re-run prepare.R with this chromosome.", chr))
gds_path <- file.path(config$gds_dir, gsub("\\{chr\\}", chr, config$gds_pattern))
stopifnot(file.exists(gds_path))

# ---- Models + (optional) STAAR context ----
pi_models <- load_PI_models(config$pi_model_dir)
b_model   <- if (!is.null(config$b_model_path)) load_B_model(config$b_model_path) else NULL
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
}

# ---- Output policy + GC calibration + SPA option ----
policy <- glow_output_policy(
  write_flat_table     = TRUE,
  write_evidence       = isTRUE(config$write_evidence),
  write_staar_detail   = isTRUE(config$write_staar_detail),
  write_result_objects = isTRUE(config$write_result_objects))
inflation_factor <- null_or(config$inflation_factor, 1.0)
z_scale <- 1 / sqrt(inflation_factor)
if (inflation_factor != 1.0)
  log_msg(sprintf("GC calibration: inflation_factor = %.4f -> z_scale = %.4f", inflation_factor, z_scale))
use_spa <- config$use_spa
if (!is.null(use_spa)) log_msg(sprintf("SPA option: use_spa = %s", use_spa))

# ---- Run the unified Stage-2 driver (region_type = "window") ----
output_format <- null_or(config$output_format, "csv")
ext <- if (identical(output_format, "fst")) "fst" else "csv"
out <- run_scan_unit(
  gds_path = gds_path, null_model = null_model, pi_models = pi_models,
  pi_features = config$pi_features, reference_medians = ref_medians,
  b_func = config$b_func, b_model = b_model,
  region_type = "window", config = config,
  results_path = file.path(results_dir, sprintf("scan_chunk_%04d.%s", chunk_id, ext)),
  output_dir = output_dir, output_policy = policy,
  chunk = chunk, unit_id = sprintf("chunk_%04d", chunk_id),
  staar_context = staar_context, use_spa = use_spa, z_scale = z_scale,
  max_regions = max_windows, progress_every = 100L, logger = log_msg, verbose = 0)

s <- out$summary
log_msg(sprintf(paste0("Done chunk %d. n_windows=%d, n_ok=%d, n_skip_empty=%d, ",
                       "n_error=%d, scan=%.1fs, wall=%.1fs"),
                chunk_id, s$n_windows, s$n_ok, s$n_skip_empty, s$n_error,
                s$runtime_secs, as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
