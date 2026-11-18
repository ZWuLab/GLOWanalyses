#!/usr/bin/env Rscript
# ============================================================================
# 01-training/pi/annotate-pi-case.R - FAVOR-annotate the PI-training CASE set
#
# Cohort-agnostic. Wraps GLOWr::annotate_favor on the case csv from
# prepare-pi-case.R, adding the FAVOR functional features that PI models train
# on, and writes the annotated case csv that train-pi.R / evaluate-pi.R consume.
#
# Where it runs: the settings come from the YAML configuration, and the output folder is
# <output_root>/<dataset>/<analysis>/training/pi/ (output_dir below). The script reads the
# case csv that prepare-pi-case.R wrote into the same folder.
#
# Outputs (under <output_dir>/):
#   pi_case_data_annotated.csv   case variants + FAVOR features (PI train input)
#
# Usage (from the GLOWanalyses directory):
#   ./glow run training --config <study>.yml --part pi_annotate     # the usual call: every unit of the part
#   Rscript 01-training/pi/annotate-pi-case.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))
suppressMessages({ library(SeqArray); library(gdsfmt) })

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Configuration: the YAML file, through the stage helper ----
source("lib/glow-stage.R")
st <- glow_stage_begin("training", "pi_annotate")
config_path <- st$config_path
cat("Config: ", config_path, "\n", sep = "")
config <- environment()
g0 <- function(nm, d = NULL) get0(nm, envir = config, ifnotfound = d, inherits = FALSE)

# ---- Run identity: the output folder and the run name, from the layout of the configuration ----
output_dir <- st$output_dir
run_name   <- st$run_name
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
na_handling  <- null_or(g0("favor_na_handling"), "keep")       # "keep" | "zero" | "drop"
db_format    <- null_or(g0("favor_db_format"), "auto")         # "auto" | "csv" | "parquet"
fav_release  <- g0("favor_release")                             # NULL -> "unknown"
rsid_pol     <- null_or(g0("favor_rsid_policy"), "require")     # "require" | "record"

# ---- Annotate (GLOWr::annotate_favor) ----
cat(sprintf("Annotating %d case variants with FAVOR [match_method = %s]...\n",
            nrow(cases), match_method))
start_time <- Sys.time()
annot_args <- list(
  variants        = cases,
  favor_db_path   = favor_db,
  match_method    = match_method,
  use_xsv         = use_xsv,
  na_handling     = na_handling,
  favor_db_format = db_format,
  favor_release   = fav_release,
  rsid_policy     = rsid_pol,
  verbose         = 2)
# favor_features = NULL means the package default set: omit the argument.
if (!is.null(features)) annot_args$features <- features
cases_annotated <- do.call(annotate_favor, annot_args)
cat(sprintf("Annotation completed in %.2f min\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))

# ---- Persist the annotated case csv ----
out_csv <- file.path(output_dir, "pi_case_data_annotated.csv")
data.table::fwrite(cases_annotated, out_csv)
cat(sprintf("Annotated case set: %d variants x %d columns -> %s\n",
            nrow(cases_annotated), ncol(cases_annotated), out_csv))
cat(sprintf("Next: train-pi.R --config %s\n", config_path))

# The completion record of this unit.
glow_stage_end(st, outputs = out_csv)
