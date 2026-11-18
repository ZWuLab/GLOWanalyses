#!/usr/bin/env Rscript
# ============================================================================
# 01-training/pi/prepare-pi-case.R - clean known-associated SNPs into a
#                                    PI-training CASE set
#
# Cohort-agnostic. Wraps GLOWr::prepare_PI_case_data: reads a raw table of known
# trait-associated SNPs, applies the config-supplied column mapping, author /
# chromosome exclusions, and QC, and writes the cleaned CASE csv that the FAVOR
# annotation step (annotate-pi-case.R) then annotates for PI model training.
#
# Where it runs: the settings come from the YAML configuration, and the output folder is
# the one of its step in the derived layout, <output_root>/<dataset>/<analysis>/training/pi/
# (output_dir below).
#
# Outputs (under <output_dir>/):
#   pi_case_data.csv   cleaned case variants (the input to annotate-pi-case.R)
#
# Usage (from the GLOWanalyses directory):
#   ./glow run training --config <study>.yml --part pi_prepare     # the usual call: every unit of the part
#   Rscript 01-training/pi/prepare-pi-case.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Configuration: the YAML file, through the stage helper ----
source("lib/glow-stage.R")
st <- glow_stage_begin("training", "pi_prepare")
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
pi_raw_path <- g0("pi_raw_path")
if (is.null(pi_raw_path) || !file.exists(pi_raw_path))
  stop("Config `pi_raw_path` must be an existing raw known-SNP table ",
       "(path to the .xlsx/.csv/.txt of known trait-associated SNPs).")

# ---- Optional config fields -> explicit defaults (forwarded to GLOWr) ----
pi_format        <- g0("pi_format", "auto")           # file-format hint for the loader
column_mapping   <- g0("pi_column_mapping")           # NULL -> loader's column auto-detect
filter_include   <- g0("pi_filter_include")           # NULL -> no include filter
filter_exclude   <- g0("pi_filter_exclude")           # NULL -> no exclude filter
exclude_authors  <- g0("pi_exclude_authors")          # NULL -> keep all authors
exclude_chroms   <- g0("pi_exclude_chromosomes")      # NULL -> keep all chromosomes
qc_filters       <- g0("pi_qc_filters")               # NULL -> prepare_PI_case_data() default QC

cat(sprintf("Raw known-SNP table: %s\n", pi_raw_path))

# ---- Prepare the case set (GLOWr::prepare_PI_case_data) ----
# Omit NULL list fields so the GLOWr default applies for each.
prep_args <- list(
  data                = pi_raw_path,
  format              = pi_format,
  column_mapping      = column_mapping,
  filter_include      = filter_include,
  filter_exclude      = filter_exclude,
  exclude_authors     = exclude_authors,
  exclude_chromosomes = exclude_chroms,
  verbose             = 1)
if (!is.null(qc_filters)) prep_args$qc_filters <- qc_filters

result <- do.call(prepare_PI_case_data, prep_args)

cat(sprintf("\nPI case set: %d variants x %d columns\n",
            nrow(result$data), ncol(result$data)))

# ---- Persist the cleaned case csv ----
out_csv <- file.path(output_dir, "pi_case_data.csv")
data.table::fwrite(result$data, out_csv)
cat("PI case set written: ", out_csv, "\n", sep = "")
cat(sprintf("Next: annotate-pi-case.R --config %s\n", config_path))

# The completion record of this unit.
glow_stage_end(st, outputs = out_csv)
