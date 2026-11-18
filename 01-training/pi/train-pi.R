#!/usr/bin/env Rscript
# ============================================================================
# 01-training/pi/train-pi.R - train the PI (prior-inclusion-probability) ensemble
#
# Cohort-agnostic. Wraps GLOWr::train_PI_models on the annotated CASE set from
# annotate-pi-case.R against a CONTROL source (the cohort's FAVOR-annotated tree,
# config-provided), writing the trained PI ensemble that the SNV-set step reads
# (`snv_set.pi_model_dir`, derived from the layout).
#
# CONTROL ANNOTATIONS are produced UPSTREAM, not here: the per-chromosome FAVOR
# aGDS/CSV control tree is built by 00-data-prep/annotate-favor.R (in-pipeline) or
# the whole-genome HPC array driver system.file("scripts/annotate_favor_hpc.sh", package = "GLOWr").
# This template only points train_PI_models at that already-annotated tree via
# `pi_control_source`.
#
# Where it runs: the settings come from the YAML configuration, and the output folder is
# <output_root>/<dataset>/<analysis>/training/pi/ (output_dir below). The script reads the
# annotated case csv from the same folder. The models are written under
# <output_dir>/<pi_models_subdir> (default "pi_models").
#
# Outputs (under <output_dir>/<pi_models_subdir>/):
#   model_*.rds              the n_models trained ensemble members
#   training_metadata.rds    train_PI_models() metadata snapshot
#
# Usage (from the GLOWanalyses directory):
#   ./glow run training --config <study>.yml --part pi_train     # the usual call: every unit of the part
#   Rscript 01-training/pi/train-pi.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))
suppressMessages(library(data.table))

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Configuration: the YAML file, through the stage helper ----
source("lib/glow-stage.R")
st <- glow_stage_begin("training", "pi_train")
config_path <- st$config_path
cat("Config: ", config_path, "\n", sep = "")
config <- environment()
g0 <- function(nm, d = NULL) get0(nm, envir = config, ifnotfound = d, inherits = FALSE)

# ---- Run identity: the output folder and the run name, from the layout of the configuration ----
output_dir <- st$output_dir
run_name   <- st$run_name
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Run: %s  ->  output_dir = %s\n", run_name, output_dir))

# ---- Annotated CASE csv (written by annotate-pi-case.R); config may override ----
case_csv <- null_or(g0("pi_case_annotated_csv"),
                    file.path(output_dir, "pi_case_data_annotated.csv"))
if (!file.exists(case_csv))
  stop("Annotated case csv not found: ", case_csv,
       "\nRun annotate-pi-case.R first (or set `pi_case_annotated_csv`).")

# ---- Required config field: CONTROL source (the cohort FAVOR-annotated tree) ----
# Built upstream by 00-data-prep/annotate-favor.R or the HPC array driver
# system.file("scripts/annotate_favor_hpc.sh", package = "GLOWr"). Pass that already-annotated
# per-chromosome aGDS/CSV directory (or a single annotation file) here.
control_source <- g0("pi_control_source")
if (is.null(control_source) || !(dir.exists(control_source) || file.exists(control_source)))
  stop("Config `pi_control_source` must be the cohort's FAVOR-annotated control ",
       "tree (a directory of per-chromosome aGDS/CSV) or a single annotation file. ",
       "Build it upstream with 00-data-prep/annotate-favor.R or the bundled GLOWr ",
       "HPC driver at system.file('scripts/annotate_favor_hpc.sh', package = 'GLOWr').")
cat("Control source (annotated upstream): ", control_source, "\n", sep = "")

# ---- Optional config fields -> explicit defaults (forwarded to train_PI_models) ----
n_models            <- as.integer(null_or(g0("pi_n_models"), 100L))
controls_multiplier <- null_or(g0("pi_controls_multiplier"), 1)
max_controls        <- g0("pi_max_controls", 50000L)   # NULL = load all controls
model_type          <- null_or(g0("pi_model_type"), "LASSO")  # "GLM" | "LASSO"
chromosomes         <- as.integer(null_or(g0("pi_chromosomes"), 1:22))
features            <- g0("pi_features")               # NULL -> .default_PI_features()
random_seed         <- as.integer(null_or(g0("pi_random_seed"), 42L))
# The format of the control files. A control folder may hold CSV and GDS files side by side,
# and "auto" would then read the CSV files. So the format that the configuration names is
# passed on, and the training reads the files that the evaluation reads and the completion
# record states. (GLOWr 0.3.6 and later; "agds" is another word for "gds".)
control_format      <- tolower(null_or(g0("pi_control_format"), "auto"))
if (control_format == "agds") control_format <- "gds"

# ---- Model output dir: an output_dir-relative subdir (the SNV-set pi_model_dir) ----
pi_models_subdir <- null_or(g0("pi_models_subdir"), "pi_models")
models_dir <- file.path(output_dir, pi_models_subdir)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)
# The SNV-set step reads every model in this folder, so the folder starts empty (a YAML
# configuration): a model that an earlier run left must not join this ensemble.
glow_stage_clear_dir(st, models_dir)

cat(sprintf(paste0("Training PI ensemble [n_models = %d, model_type = %s, ",
                   "controls_multiplier = %s, max_controls = %s, seed = %d]\n"),
            n_models, model_type, format(controls_multiplier),
            if (is.null(max_controls)) "all" else format(max_controls), random_seed))

# ---- Train (GLOWr::train_PI_models) ----
result <- train_PI_models(
  case_csv            = case_csv,
  control_source      = control_source,
  output_dir          = models_dir,
  n_models            = n_models,
  controls_multiplier = controls_multiplier,
  max_controls        = max_controls,
  model_type          = model_type,
  features            = features,
  chromosomes         = chromosomes,
  random_seed         = random_seed,
  verbose             = 1,
  control_format      = control_format)

# ---- Persist the metadata snapshot + report ----
metadata_file <- file.path(models_dir, "training_metadata.rds")
saveRDS(result$metadata, metadata_file)
n_rds <- length(list.files(models_dir, pattern = "\\.rds$"))
cat(sprintf("PI ensemble trained: %d model(s) -> %s\n",
            length(result$models), models_dir))
cat(sprintf("  metadata -> %s ; total .rds files in dir = %d\n", metadata_file, n_rds))
cat(sprintf("Next: evaluate-pi.R --config %s\n", config_path))

# The completion record of this unit.
glow_stage_end(st, outputs = models_dir)
