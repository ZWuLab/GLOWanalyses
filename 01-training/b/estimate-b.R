#!/usr/bin/env Rscript
# ============================================================================
# 01-training/b/estimate-b.R - fit a B (effect-size distribution) model
#
# Cohort-agnostic. Wraps GLOWr::train_B_model on the cleaned training set from
# prepare-b-training-data.R, persists it with GLOWr::save_B_model, and (when the
# config supplies a target MAF vector) evaluates B(MAF) with GLOWr::predict_B
# (equivalently GLOWr::get_B) for inspection. The trained b_model.rds is the
# B source that the SNV-set step reads (`snv_set.b_model_path`, derived from the layout).
#
# Where it runs: the settings come from the YAML configuration, and the output folder is
# <output_root>/<dataset>/<analysis>/training/b/ (output_dir below). The script reads the
# training set that prepare-b-training-data.R wrote into the same folder.
#
# Outputs (under <output_dir>/):
#   b_model.rds   trained B model (save_B_model format; the SNV-set b_model_path)
#   b_values.csv  (only if `b_target_maf` set) target MAF -> predicted B table
#
# Usage (from the GLOWanalyses directory):
#   ./glow run training --config <study>.yml --part b_estimate     # the usual call: every unit of the part
#   Rscript 01-training/b/estimate-b.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Configuration: the YAML file, through the stage helper ----
source("lib/glow-stage.R")
st <- glow_stage_begin("training", "b_estimate")
config_path <- st$config_path
cat("Config: ", config_path, "\n", sep = "")
config <- environment()
g0 <- function(nm, d = NULL) get0(nm, envir = config, ifnotfound = d, inherits = FALSE)

# ---- Run identity: the output folder and the run name, from the layout of the configuration ----
output_dir <- st$output_dir
run_name   <- st$run_name
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
  # A null b_target_trait means the training trait, as the fitted model holds it: a null
  # b_training_trait is "mixed" there (GLOWr::train_B_model()). (Through release 0.3.6 a null
  # was passed on as null, because the field exists in the base configuration with the value
  # NULL.)
  training_trait_used <- null_or(b_model$training_summary$trait_type, null_or(training_trait, "mixed"))
  target_trait     <- null_or(g0("b_target_trait"), training_trait_used)
  target_case_prop <- g0("b_target_case_prop")  # binary-only; NULL for continuous
  # The trait type matters for a model fitted by the p-value method only: GLOWr::predict_B()
  # converts its prediction (the square root of a heritability) to B by the trait type. The
  # pairs that cannot be converted stop here, with the settings named. A target MAF of zero
  # has B = 0 whatever the trait, so targets that are all zero need no conversion.
  b_method_used <- b_model$method_used
  if (identical(b_method_used, "both")) b_method_used <- b_model$comparison$method_selected
  if (identical(b_method_used, "pvalue_method") && any(target_maf != 0)) {
    if (identical(target_trait, "mixed"))
      stop("The B table of `b_target_maf` needs a binary or continuous trait: the B model was fitted ",
           "by the p-value method, and the trait is \"mixed\" (`b_target_trait`, or `b_training_trait` ",
           "when that is null). Set `b_target_trait` to binary or continuous.")
    if (identical(target_trait, "binary") && is.null(target_case_prop))
      stop("The B table of `b_target_maf` for a binary trait needs `b_target_case_prop`, the ",
           "proportion of cases in the target study: the B model was fitted by the p-value method.")
    if (identical(target_trait, "continuous"))
      stop("The B table of `b_target_maf` for a continuous trait and a B model fitted by the p-value ",
           "method needs the standard error of the trait, which the configuration does not hold. ",
           "Set `b_target_maf` to null and call GLOWr::predict_B(target_SE = ) on the saved model.")
  }
  B_pred <- predict_B(b_model, target_maf,
                      target_trait     = if (identical(target_trait, "mixed")) NULL else target_trait,
                      target_case_prop = target_case_prop)
  out_b_values <- file.path(output_dir, "b_values.csv")
  write.csv(data.frame(target_MAF = target_maf, B = B_pred),
            out_b_values, row.names = FALSE)
  cat(sprintf("Predicted B for %d target MAF(s) [range %.4g, %.4g] -> %s\n",
              length(target_maf), min(B_pred), max(B_pred), out_b_values))
}

cat(sprintf("B estimation complete for run '%s'.\n", run_name))

# The completion record of this unit.
# The table of B values is an output when `b_target_maf` asked for it. Without it, a table that
# an earlier run left is removed, so that the folder holds the files of this run only.
b_values_csv <- file.path(output_dir, "b_values.csv")
if (is.null(target_spec)) unlink(b_values_csv, expand = FALSE)
glow_stage_end(st, outputs = c(out_model, if (!is.null(target_spec)) b_values_csv))
