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
# RUN-ORG (mirrors 03-snv-set): output_dir is DERIVED from the --config path;
# run_name = that dir's name. The run config SOURCES base-config/training_base.R.
#
# Outputs (under <output_dir>/):
#   pi_case_data.csv   cleaned case variants (the input to annotate-pi-case.R)
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 01-training/pi/prepare-pi-case.R \
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
