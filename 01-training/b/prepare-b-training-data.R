#!/usr/bin/env Rscript
# ============================================================================
# 01-training/b/prepare-b-training-data.R - clean a known-effect SNP table
#                                           into a B-training set
#
# Cohort-agnostic. Wraps GLOWr::prepare_B_training_data: reads a raw table of
# known trait-associated SNPs (effect sizes / p-values / MAF / N), applies the
# config-supplied column mapping, exclusion/inclusion filters, and QC, and
# writes the cleaned training set (the input every B model fit consumes) to the
# run's outputs/.
#
# RUN-ORG (mirrors 03-snv-set): a RUN is a directory holding config.R and a
# sibling outputs/. output_dir is DERIVED from the --config path (never declared
# in a config), and run_name = that dir's name. The run config SOURCES
# base-config/training_base.R, then sets the raw known-SNP path + any overrides.
#
# Outputs (under <output_dir>/):
#   b_training_data.rds   full prepare_B_training_data() result (data + metadata)
#   b_training_data.csv   the cleaned training table only (for inspection)
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 01-training/b/prepare-b-training-data.R \
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
