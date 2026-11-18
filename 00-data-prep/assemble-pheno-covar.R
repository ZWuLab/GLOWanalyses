#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/assemble-pheno-covar.R - build the (sample_id, Y, X, trait) bundle
#
# Cohort-agnostic. Wraps GLOWr::assemble_pheno_covar: reads the sample order
# from one chromosome's GDS, an external phenotype/covariate table, and a PC
# table, and writes the analysis bundle that fit_null_model / marginal_scan /
# the SNV-set pipeline (prepare_scan_run) consume, into the data_prep folder of the analysis,
# <output_root>/<dataset>/<analysis>/data_prep/. One bundle per analysis (its covariates).
#
# Usage (from the GLOWanalyses directory):
#   ./glow run data_prep --config <study>.yml --part assemble_pheno     # the usual call: every unit of the part
#   Rscript 00-data-prep/assemble-pheno-covar.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))
suppressMessages(library(SeqArray))
source("00-data-prep/_dataprep_lib.R")
source("lib/glow-stage.R")
st  <- glow_stage_begin("data_prep", "assemble_pheno")
pa  <- list(config = st$config_path, opts = st$opts)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)
run_name <- st$run_name

match_meth <- null_or(g0("favor_match_method"), "flexible")
# The folders come from the layout (stage_cohort_paths()), never rebuilt here.
paths      <- stage_cohort_paths(st, g0, match_method = match_meth)

# ---- Sample order from one chromosome's (a)GDS ----
sample_chr <- as.integer(null_or(g0("sample_gds_chr"), 22L))
# The annotated GDS files are the user's when the configuration names their folder
# (shared.analysis_gds_dir with shared.analysis_gds_pattern), and those of the run otherwise.
favor_pat  <- null_or(paths$read_favor_gds_pattern, null_or(g0("favor_gds_pattern"), "chr{chr}_hg38_favor.gds"))
sample_gds <- chr_path(null_or(paths$read_favor_gds_dir, paths$favor_gds_dir), favor_pat, sample_chr)
if (!file.exists(sample_gds))
  stop("The GDS file that gives the sample order does not exist: ", sample_gds,
       " (chromosome `sample_gds_chr` of the annotated GDS files).")
gds <- seqOpen(sample_gds, readonly = TRUE)
sample_ids <- seqGetData(gds, "sample.id")

# ---- Outcome values (from a GDS sample.annotation node OR a pheno column) ----
outcome_cfg <- g0("outcome")
if (is.null(outcome_cfg) || (is.null(outcome_cfg$node) && is.null(outcome_cfg$col)))
  stop("Config `outcome` must give a `node` (GDS sample.annotation) or `col` (pheno column), plus `map`.")
pheno_csv  <- g0("pheno_csv"); pheno_id_col <- g0("pheno_id_col")
stopifnot(!is.null(pheno_csv), file.exists(pheno_csv), !is.null(pheno_id_col))
pheno <- as.data.frame(data.table::fread(pheno_csv))

if (!is.null(outcome_cfg$node)) {
  outcome_values <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(
    gds, paste0("sample.annotation/", outcome_cfg$node)))      # already in GDS sample order
} else {
  oidx <- match(sample_ids, as.character(pheno[[pheno_id_col]]))
  outcome_values <- pheno[[outcome_cfg$col]][oidx]
}
seqClose(gds)

# ---- PC table ----
pcs_cfg <- g0("pcs")
pcs_arg <- NULL
if (!is.null(pcs_cfg)) {
  stopifnot(!is.null(pcs_cfg$path), file.exists(pcs_cfg$path))
  scores <- if (grepl("\\.rds$", pcs_cfg$path, ignore.case = TRUE)) readRDS(pcs_cfg$path)
            else as.data.frame(data.table::fread(pcs_cfg$path))
  # `cols` null takes the first `n_pcs` components: the first columns of the table that
  # are not the identifier column.
  pcs_cols <- pcs_cfg$cols
  if (is.null(pcs_cols)) {
    n_use <- as.integer(null_or(g0("n_pcs"), 20L))
    candidates <- setdiff(colnames(scores), pcs_cfg$id_col)
    if (length(candidates) < n_use)
      stop("`pcs$cols` is null, which takes the first `n_pcs` = ", n_use, " components, and the PC table has ",
           length(candidates), " column(s) besides the identifier: ", pcs_cfg$path)
    pcs_cols <- candidates[seq_len(n_use)]
  }
  pcs_arg <- list(scores = scores, id_col = pcs_cfg$id_col, cols = pcs_cols)
}

# ---- Assemble ----
bundle <- assemble_pheno_covar(
  sample_ids   = sample_ids,
  outcome      = list(values = outcome_values, map = outcome_cfg$map),
  trait        = null_or(g0("trait"), "binary"),
  pheno        = pheno, pheno_id_col = pheno_id_col,
  covariates   = g0("covariates"),
  pcs          = pcs_arg,
  drop_incomplete = isTRUE(null_or(g0("drop_incomplete"), TRUE)),
  covar_set    = run_name, verbose = 1)

dir.create(paths$pheno_dir, recursive = TRUE, showWarnings = FALSE)
# The bundle: <analysis>/data_prep/<pheno_out_name>, with links to the dataset's data folders
# beside it.
pheno_out_name <- null_or(g0("pheno_out_name"), "pheno_covar.rds")
out_rds <- file.path(paths$pheno_dir, pheno_out_name)
glow_link_dataset_assets(st$layout)
saveRDS(bundle, out_rds)
write_dataprep_provenance(
  paths$pheno_dir, step = "assemble_pheno_covar",
  script = "00-data-prep/assemble-pheno-covar.R",
  source_alias = basename(pheno_csv), source_path = pheno_csv,
  config_snapshot = list(run = run_name, outcome = outcome_cfg, covariates = g0("covariates"),
                         pcs = pcs_cfg, trait = g0("trait")))
cat(sprintf("Bundle: %d samples, %d covariates [%s] -> %s\n",
            length(bundle$Y), ncol(bundle$X), paste(bundle$covar_names, collapse = ", "), out_rds))

# The completion record of this unit.
glow_stage_end(st, outputs = out_rds)
