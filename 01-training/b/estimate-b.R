#!/usr/bin/env Rscript
# ============================================================================
# 01-training/b/estimate-b.R - fit a B (effect-size distribution) model
#
# Cohort-agnostic. Wraps GLOWr::train_B_model on the cleaned training set from
# prepare-b-training-data.R, persists it with GLOWr::save_B_model, and (when the
# config supplies a target MAF vector) evaluates B(MAF) with GLOWr::predict_B
# (equivalently GLOWr::get_B) for inspection. The trained b_model.rds is the
# B source an SNV-set run points at via `b_model_path` (snv_set_base.R).
#
# RUN-ORG (mirrors 03-snv-set): output_dir is DERIVED from the --config path;
# this script reads the training set written by prepare-b-training-data.R in the
# SAME run's outputs/. The run config SOURCES base-config/training_base.R.
#
# Outputs (under <output_dir>/):
#   b_model.rds   trained B model (save_B_model format; the SNV-set b_model_path)
#   b_values.csv  (only if `b_target_maf` set) target MAF -> predicted B table
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 01-training/b/estimate-b.R \
#       --config runs/<name>/config.R

suppressMessages(library(GLOWr))

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <run>/config.R (required) ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) {
    config_path <- args[i + 1L]; i <- i + 2L
  } else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))
cat("Config: ", config_path, "\n", sep = "")

# ---- Source config (VALUES only) ----
source(config_path, local = TRUE)
config <- environment()
g0 <- function(nm, d = NULL) get0(nm, envir = config, ifnotfound = d, inherits = FALSE)

# ---- Run identity: output_dir + run_name DERIVED from the --config location ----
if (!is.null(g0("output_dir")))
  message("Note: config set `output_dir`; ignoring it (output_dir is derived ",
          "from the --config path, not declared).")
output_dir <- file.path(dirname(config_path), "outputs")
run_name   <- basename(dirname(config_path))
# Shared-output OPT-IN: NULL (default) keeps outputs/ local; a shared-root path
# born-shares this run's outputs/ (symlink into <root>/<repo-relative-path>) before
# it is populated. GLOWpipeline is loaded only on opt-in, so default local runs stay
ssr <- get0("symlinked_shared_root", envir = config, ifnotfound = NULL, inherits = FALSE)
if (!is.null(ssr)) {
  suppressMessages(library(GLOWpipeline))
  GLOWpipeline::ensure_shared_output_dir(output_dir, share_root = ssr)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Run: %s  ->  output_dir = %s\n", run_name, output_dir))

# ---- Load the cleaned training set (written by prepare-b-training-data.R) ----
# A config may override the location (b_training_rds); default = this run's output.
b_training_rds <- null_or(g0("b_training_rds"), file.path(output_dir, "b_training_data.rds"))
if (!file.exists(b_training_rds))
  stop("B training set not found: ", b_training_rds,
       "\nRun prepare-b-training-data.R first (or set `b_training_rds`).")
training <- readRDS(b_training_rds)$data
cat(sprintf("Training set: %d variants from %s\n", nrow(training), b_training_rds))

# ---- Config: training trait, method, and which training columns to use ----
training_trait <- g0("b_training_trait", "binary")  # "binary" | "continuous" | "mixed"
method         <- g0("b_method", "auto")            # "beta" | "pvalue" | "both" | "auto"

# Column names IN the training table (default to the prepare_B_training_data
# output schema MAF/BETA/P/N). NULL fields are dropped so train_B_model picks the
# method that the available columns support.
maf_col  <- null_or(g0("b_maf_col"),  "MAF")
beta_col <- g0("b_beta_col", "BETA")  # NULL skips BETA (no beta method)
p_col    <- g0("b_p_col",    "P")     # NULL skips P    (no pvalue method)
n_col    <- g0("b_n_col",    "N")     # NULL skips N    (no pvalue method)

if (!maf_col %in% names(training))
  stop("MAF column '", maf_col, "' not in the training set (have: ",
       paste(names(training), collapse = ", "), ").")

train_args <- list(
  training_trait = training_trait,
  training_MAF   = training[[maf_col]],
  method         = method,
  verbose        = 1)
if (!is.null(beta_col)) {
  if (!beta_col %in% names(training)) stop("BETA column '", beta_col, "' not in the training set.")
  train_args$training_BETA <- training[[beta_col]]
}
if (!is.null(p_col)) {
  if (!p_col %in% names(training)) stop("P column '", p_col, "' not in the training set.")
  train_args$training_P <- training[[p_col]]
}
if (!is.null(n_col)) {
  if (!n_col %in% names(training)) stop("N column '", n_col, "' not in the training set.")
  train_args$training_N <- training[[n_col]]
}

cat(sprintf("Fitting B model [trait = %s, method = %s] on columns: %s\n",
            training_trait, method,
            paste(intersect(c("MAF", "BETA", "P", "N"),
                            c(if (!is.null(maf_col)) "MAF",
                              if (!is.null(beta_col)) "BETA",
                              if (!is.null(p_col)) "P",
                              if (!is.null(n_col)) "N")), collapse = ", ")))

# ---- Fit + persist the B model (GLOWr::train_B_model -> save_B_model) ----
b_model  <- do.call(train_B_model, train_args)
out_model <- file.path(output_dir, "b_model.rds")
save_B_model(b_model, out_model)
cat("B model written: ", out_model, "\n", sep = "")

# ---- Optional: evaluate B for a config-supplied target MAF vector ----
# `b_target_maf` may be a numeric vector or a length-3 spec list(from, to, n) for
# an evenly-spaced grid. When set we write b_values.csv via predict_B.
target_spec <- g0("b_target_maf")
if (!is.null(target_spec)) {
  if (is.list(target_spec)) {
    stopifnot(!is.null(target_spec$from), !is.null(target_spec$to))
    n_grid     <- as.integer(null_or(target_spec$n, 500L))
    target_maf <- seq(target_spec$from, target_spec$to, length.out = n_grid)
  } else {
    target_maf <- as.numeric(target_spec)
  }
  target_trait     <- g0("b_target_trait", training_trait)
  target_case_prop <- g0("b_target_case_prop")  # binary-only; NULL for continuous
  B_pred <- predict_B(b_model, target_maf,
                      target_trait     = target_trait,
                      target_case_prop = target_case_prop)
  out_b_values <- file.path(output_dir, "b_values.csv")
  write.csv(data.frame(target_MAF = target_maf, B = B_pred),
            out_b_values, row.names = FALSE)
  cat(sprintf("Predicted B for %d target MAF(s) [range %.4g, %.4g] -> %s\n",
              length(target_maf), min(B_pred), max(B_pred), out_b_values))
}

cat(sprintf("B estimation complete for run '%s'.\n", run_name))
