#!/usr/bin/env Rscript
# ============================================================================
# 01-training/pi/evaluate-pi.R - evaluate a trained PI ensemble (AUC + ROC)
#
# Cohort-agnostic. Loads the run's trained PI ensemble + the annotated CASE set
# and a CONTROL sample (the cohort FAVOR-annotated tree, config-provided), then
# wraps GLOWr::evaluate_PI_models + GLOWr::plot_PI_roc to write a per-model AUC
# table and a ROC pdf. Mirrors the core of paper1-16-pi-training/05 for ONE model set.
#
# Where it runs: the settings come from the YAML configuration, and the output folder is
# <output_root>/<dataset>/<analysis>/training/pi/ (output_dir below). The script reads the
# annotated case csv and the trained ensemble from the same folder.
#
# Outputs (under <output_dir>/):
#   pi_eval_auc_per_model.csv    per-model AUC for the ensemble
#   pi_eval_summary.csv          mean/sd/min/max/ensemble AUC (one row)
#   pi_eval_roc.pdf              ROC curves (individual + ensemble)
#
# Usage (from the GLOWanalyses directory):
#   ./glow run training --config <study>.yml --part pi_evaluate     # the usual call: every unit of the part
#   Rscript 01-training/pi/evaluate-pi.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))
suppressMessages({ library(data.table); library(Matrix); library(glmnet) })

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Configuration: the YAML file, through the stage helper ----
source("lib/glow-stage.R")
st <- glow_stage_begin("training", "pi_evaluate")
config_path <- st$config_path
cat("Config: ", config_path, "\n", sep = "")
config <- environment()
g0 <- function(nm, d = NULL) get0(nm, envir = config, ifnotfound = d, inherits = FALSE)

# ---- Run identity: the output folder and the run name, from the layout of the configuration ----
output_dir <- st$output_dir
run_name   <- st$run_name
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Run: %s  ->  output_dir = %s\n", run_name, output_dir))

# ---- Inputs from the SAME run (config may override each path) ----
case_csv <- null_or(g0("pi_case_annotated_csv"),
                    file.path(output_dir, "pi_case_data_annotated.csv"))
if (!file.exists(case_csv))
  stop("Annotated case csv not found: ", case_csv,
       "\nRun annotate-pi-case.R first (or set `pi_case_annotated_csv`).")

pi_models_subdir <- null_or(g0("pi_models_subdir"), "pi_models")
models_dir <- null_or(g0("pi_models_dir"), file.path(output_dir, pi_models_subdir))
if (!dir.exists(models_dir))
  stop("Trained PI ensemble dir not found: ", models_dir,
       "\nRun train-pi.R first (or set `pi_models_dir`).")

control_source <- g0("pi_control_source")
if (is.null(control_source) || !(dir.exists(control_source) || file.exists(control_source)))
  stop("Config `pi_control_source` must be the cohort's FAVOR-annotated control ",
       "tree (a directory of per-chromosome aGDS/CSV) or a single annotation file.")

# ---- Optional config fields -> explicit defaults ----
features          <- null_or(g0("pi_features"), .default_PI_features())
data_type         <- null_or(g0("pi_control_format"), "auto")   # "auto" | "gds" | "csv"
eval_max_controls <- g0("pi_eval_max_controls", 5000L)          # NULL = all controls
n_eval_controls   <- g0("pi_n_eval_controls")                   # NULL = use max_controls
chromosomes       <- as.integer(null_or(g0("pi_chromosomes"), 1:22))
random_seed       <- as.integer(null_or(g0("pi_random_seed"), 42L))

# ---- Load the trained ensemble ----
cat("Loading PI ensemble from ", models_dir, " ...\n", sep = "")
model_result <- load_PI_models(models_dir)
cat(sprintf("  %d %s model(s) loaded\n", model_result$n_models, model_result$model_type))

# ---- Load case annotations (matrix of the evaluated features) ----
cases <- fread(case_csv, data.table = FALSE)
missing_feats <- setdiff(features, names(cases))
if (length(missing_feats) > 0)
  stop("Case csv missing PI feature(s): ", paste(missing_feats, collapse = ", "))
case_mat <- as.matrix(cases[, features, drop = FALSE])
cat(sprintf("  Cases: %d variants x %d features\n", nrow(case_mat), ncol(case_mat)))

# ---- Load a control sample (FAVOR-annotated tree built upstream) ----
cat("Loading control annotations [", data_type, "] from ", control_source, " ...\n", sep = "")
set.seed(random_seed)
ctrl_mat <- load_control_annotations(
  source       = control_source,
  features     = features,
  chromosomes  = chromosomes,
  format       = data_type,
  max_controls = eval_max_controls)
cat(sprintf("  Controls: %d variants x %d features\n", nrow(ctrl_mat), ncol(ctrl_mat)))

# ---- Evaluate (GLOWr::evaluate_PI_models) ----
eval_result <- evaluate_PI_models(
  models          = model_result$models,
  cases           = case_mat,
  controls        = ctrl_mat,
  model_type      = model_result$model_type,
  features        = features,
  max_controls    = NULL,                # already limited during loading
  n_eval_controls = n_eval_controls,
  random_seed     = random_seed)

cat(sprintf("Mean individual AUC = %.4f | Ensemble AUC = %.4f\n",
            eval_result$summary$mean_auc, eval_result$ensemble$auc))

# ---- Write the AUC tables ----
auc_per_model_csv <- file.path(output_dir, "pi_eval_auc_per_model.csv")
fwrite(eval_result$per_model, auc_per_model_csv)

summary_csv <- file.path(output_dir, "pi_eval_summary.csv")
fwrite(data.frame(
  run          = run_name,
  n_models     = model_result$n_models,
  model_type   = model_result$model_type,
  n_cases      = eval_result$metadata$n_cases,
  n_controls   = eval_result$metadata$n_controls,
  mean_auc     = round(eval_result$summary$mean_auc, 4),
  sd_auc       = round(eval_result$summary$sd_auc, 4),
  min_auc      = round(eval_result$summary$min_auc, 4),
  max_auc      = round(eval_result$summary$max_auc, 4),
  ensemble_auc = round(eval_result$ensemble$auc, 4),
  stringsAsFactors = FALSE), summary_csv)

# ---- ROC pdf (GLOWr::plot_PI_roc) ----
roc_pdf <- file.path(output_dir, "pi_eval_roc.pdf")
pdf(roc_pdf, width = 6, height = 6)
plot_PI_roc(eval_result, show_individual = TRUE, highlight_ensemble = TRUE,
            main = sprintf("%s PI ensemble\n(Ensemble AUC = %.3f)",
                           run_name, eval_result$ensemble$auc))
dev.off()

cat("PI evaluation written:\n",
    "  ", auc_per_model_csv, "\n",
    "  ", summary_csv, "\n",
    "  ", roc_pdf, "\n", sep = "")
cat(sprintf("PI evaluation complete for run '%s'.\n", run_name))

# The completion record of this unit.
glow_stage_end(st, outputs = c(auc_per_model_csv, summary_csv, roc_pdf))
