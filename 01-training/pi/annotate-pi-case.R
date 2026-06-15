#!/usr/bin/env Rscript
# ============================================================================
# 01-training/pi/annotate-pi-case.R - FAVOR-annotate the PI-training CASE set
#
# Cohort-agnostic. Wraps GLOWr::annotate_favor on the case csv from
# prepare-pi-case.R, adding the FAVOR functional features that PI models train
# on, and writes the annotated case csv that train-pi.R / evaluate-pi.R consume.
#
# RUN-ORG (mirrors 03-snv-set): output_dir is DERIVED from the --config path;
# this script reads the case csv written by prepare-pi-case.R in the SAME run's
# outputs/. The run config SOURCES base-config/training_base.R.
#
# Outputs (under <output_dir>/):
#   pi_case_data_annotated.csv   case variants + FAVOR features (PI train input)
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 01-training/pi/annotate-pi-case.R \
#       --config runs/<name>/config.R

suppressMessages(library(GLOWr))
suppressMessages({ library(SeqArray); library(gdsfmt) })

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

# ---- Required config field: FAVOR database ----
favor_db <- g0("favor_db")
if (is.null(favor_db) || !dir.exists(favor_db))
  stop("Config `favor_db` must be an existing FAVOR DB dir ",
       "(e.g. .../FAVOR_annotation/Essential_database_hg38).")

# ---- Input case csv (written by prepare-pi-case.R); config may override path ----
pi_case_csv <- null_or(g0("pi_case_csv"), file.path(output_dir, "pi_case_data.csv"))
if (!file.exists(pi_case_csv))
  stop("PI case csv not found: ", pi_case_csv,
       "\nRun prepare-pi-case.R first (or set `pi_case_csv`).")
cases <- data.table::fread(pi_case_csv)
cat(sprintf("Loaded %d case variants from %s\n", nrow(cases), pi_case_csv))

# ---- Optional config fields -> explicit defaults ----
match_method <- null_or(g0("favor_match_method"), "flexible")  # "exact" | "flexible"
features     <- g0("favor_features")                           # NULL -> annotate_favor() default set
use_xsv      <- isTRUE(null_or(g0("favor_use_xsv"), TRUE))
na_handling  <- null_or(g0("favor_na_handling"), "keep")       # "keep" | "remove" | "impute"

# ---- Annotate (GLOWr::annotate_favor) ----
cat(sprintf("Annotating %d case variants with FAVOR [match_method = %s]...\n",
            nrow(cases), match_method))
start_time <- Sys.time()
cases_annotated <- annotate_favor(
  variants      = cases,
  favor_db_path = favor_db,
  features      = if (is.null(features)) eval(formals(annotate_favor)$features) else features,
  match_method  = match_method,
  use_xsv       = use_xsv,
  na_handling   = na_handling,
  verbose       = 2)
cat(sprintf("Annotation completed in %.2f min\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))

# ---- Persist the annotated case csv ----
out_csv <- file.path(output_dir, "pi_case_data_annotated.csv")
data.table::fwrite(cases_annotated, out_csv)
cat(sprintf("Annotated case set: %d variants x %d columns -> %s\n",
            nrow(cases_annotated), ncol(cases_annotated), out_csv))
cat(sprintf("Next: train-pi.R --config %s\n", config_path))
