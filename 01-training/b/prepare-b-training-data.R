#!/usr/bin/env Rscript
# ============================================================================
# 01-training/b/prepare-b-training-data.R - clean a known-effect SNP table
#                                           into a B-training set
#
# Cohort-agnostic. Wraps GLOWr::prepare_B_training_data: reads a raw table of
# known trait-associated SNPs (effect sizes / p-values / MAF / N), applies the
# config-supplied column mapping, exclusion/inclusion filters, and QC, and
# writes the cleaned training set (the input every B model fit consumes) to the
# training/b folder of the analysis.
#
# Where it runs: the settings come from the YAML configuration (the raw table of known SNPs
# is `training.b_raw_path`), and the output folder is the one of its step in the derived
# layout, <output_root>/<dataset>/<analysis>/training/b/ (output_dir below).
#
# Outputs (under <output_dir>/):
#   b_training_data.rds   full prepare_B_training_data() result (data + metadata)
#   b_training_data.csv   the cleaned training table only (for inspection)
#
# Usage (from the GLOWanalyses directory):
#   ./glow run training --config <study>.yml --part b_prepare     # the usual call: every unit of the part
#   Rscript 01-training/b/prepare-b-training-data.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Configuration: the YAML file, through the stage helper ----
source("lib/glow-stage.R")
st <- glow_stage_begin("training", "b_prepare")
config_path <- st$config_path
cat("Config: ", config_path, "\n", sep = "")
config <- environment()
g0 <- function(nm, d = NULL) get0(nm, envir = config, ifnotfound = d, inherits = FALSE)

# ---- Run identity: the output folder and the run name, from the layout of the configuration ----
output_dir <- st$output_dir
run_name   <- st$run_name
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Run: %s  ->  output_dir = %s\n", run_name, output_dir))

# ---- Required config field: raw known-SNP table ----
b_raw_path <- g0("b_raw_path")
if (is.null(b_raw_path) || !file.exists(b_raw_path))
  stop("Config `b_raw_path` must be an existing raw known-SNP table ",
       "(path to the .xlsx/.csv/.txt of known trait-associated SNPs).")

# ---- Optional config fields -> explicit defaults (forwarded to GLOWr) ----
trait_type     <- g0("b_trait_type", "binary")    # "binary" | "continuous" | "mixed"
b_format       <- g0("b_format", "auto")          # file-format hint for the loader
column_mapping <- g0("b_column_mapping")          # NULL -> loader's column auto-detect
filter_include <- g0("b_filter_include")          # NULL -> no include filter
filter_exclude <- g0("b_filter_exclude")          # NULL -> no exclude filter
filter_pattern <- isTRUE(null_or(g0("b_filter_pattern"), TRUE))  # substring vs exact match
filter_autosomes <- isTRUE(null_or(g0("b_filter_autosomes"), TRUE))
qc_filters     <- g0("b_qc_filters")              # NULL -> prepare_B_training_data() default QC

cat(sprintf("Raw known-SNP table: %s [trait_type = %s]\n", b_raw_path, trait_type))

# ---- Prepare the B training set (GLOWr::prepare_B_training_data) ----
# When an optional list field is NULL we omit it so the GLOWr default applies.
prep_args <- list(
  data             = b_raw_path,
  format           = b_format,
  trait_type       = trait_type,
  column_mapping   = column_mapping,
  filter_include   = filter_include,
  filter_exclude   = filter_exclude,
  filter_pattern   = filter_pattern,
  filter_autosomes = filter_autosomes,
  verbose          = 2)
if (!is.null(qc_filters)) prep_args$qc_filters <- qc_filters

b_training <- do.call(prepare_B_training_data, prep_args)

cat(sprintf("\nB training set: %d variants x %d columns [%s]\n",
            nrow(b_training$data), ncol(b_training$data),
            paste(names(b_training$data), collapse = ", ")))

# ---- Persist: full result (.rds) + the cleaned table (.csv) ----
out_rds <- file.path(output_dir, "b_training_data.rds")
out_csv <- file.path(output_dir, "b_training_data.csv")
saveRDS(b_training, out_rds)
data.table::fwrite(b_training$data, out_csv)
cat("B training set written:\n  ", out_rds, "\n  ", out_csv, "\n", sep = "")
cat(sprintf("Next: estimate-b.R --config %s\n", config_path))

# The completion record of this unit.
glow_stage_end(st, outputs = c(out_rds, out_csv))
