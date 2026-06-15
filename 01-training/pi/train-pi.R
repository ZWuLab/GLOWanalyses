#!/usr/bin/env Rscript
# ============================================================================
# 01-training/pi/train-pi.R - train the PI (prior-inclusion-probability) ensemble
#
# Cohort-agnostic. Wraps GLOWr::train_PI_models on the annotated CASE set from
# annotate-pi-case.R against a CONTROL source (the cohort's FAVOR-annotated tree,
# config-provided), writing the trained PI ensemble that an SNV-set run points at
# via `pi_model_dir` (snv_set_base.R).
#
# CONTROL ANNOTATIONS are produced UPSTREAM, not here: the per-chromosome FAVOR
# aGDS/CSV control tree is built by 00-data-prep/annotate-favor.R (in-pipeline) or
# the whole-genome HPC array driver system.file("scripts/annotate_favor_hpc.sh", package = "GLOWr").
# This template only points train_PI_models at that already-annotated tree via
# `pi_control_source`.
#
# RUN-ORG (mirrors 03-snv-set): output_dir is DERIVED from the --config path;
# this script reads the annotated case csv from the SAME run's outputs/. Models
# are written under <output_dir>/<pi_models_subdir> (default "pi_models"). The run
# config SOURCES base-config/training_base.R.
#
# Outputs (under <output_dir>/<pi_models_subdir>/):
#   model_*.rds              the n_models trained ensemble members
#   training_metadata.rds    train_PI_models() metadata snapshot
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 01-training/pi/train-pi.R \
#       --config runs/<name>/config.R

suppressMessages(library(GLOWr))
suppressMessages(library(data.table))

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

# ---- Model output dir: an output_dir-relative subdir (the SNV-set pi_model_dir) ----
pi_models_subdir <- null_or(g0("pi_models_subdir"), "pi_models")
models_dir <- file.path(output_dir, pi_models_subdir)
dir.create(models_dir, recursive = TRUE, showWarnings = FALSE)

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
  verbose             = 1)

# ---- Persist the metadata snapshot + report ----
metadata_file <- file.path(models_dir, "training_metadata.rds")
saveRDS(result$metadata, metadata_file)
n_rds <- length(list.files(models_dir, pattern = "\\.rds$"))
cat(sprintf("PI ensemble trained: %d model(s) -> %s\n",
            length(result$models), models_dir))
cat(sprintf("  metadata -> %s ; total .rds files in dir = %d\n", metadata_file, n_rds))
cat(sprintf("Next: evaluate-pi.R --config %s\n", config_path))
