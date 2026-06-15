#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/assemble-pheno-covar.R - build the (sample_id, Y, X, trait) bundle
#
# Cohort-agnostic. Wraps GLOWr::assemble_pheno_covar: reads the sample order
# from one chromosome's GDS, an external phenotype/covariate table, and a PC
# table, and writes the analysis bundle that fit_null_model / marginal_scan /
# the SNV-set pipeline (prepare_scan_run) consume, into <base_name>_pheno/.
# One bundle per run (covariate config) -> <run>_pheno_covar.rds.
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 00-data-prep/assemble-pheno-covar.R --config <run>/config.R

suppressMessages(library(GLOWr))
suppressMessages(library(SeqArray))
source("00-data-prep/_dataprep_lib.R")

pa  <- parse_dataprep_args()
source(pa$config, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)
run_name <- basename(dirname(pa$config))

data_root <- g0("data_root"); base_name <- g0("base_name")
if (is.null(data_root) || is.null(base_name))
  stop("Config must set `data_root` and `base_name`.")
match_meth <- null_or(g0("favor_match_method"), "flexible")
paths      <- resolve_cohort_paths(data_root, base_name)

# ---- Sample order from one chromosome's (a)GDS ----
sample_chr <- as.integer(null_or(g0("sample_gds_chr"), 22L))
favor_pat  <- null_or(g0("favor_gds_pattern"), "chr{chr}_hg38_favor.gds")
sample_gds <- chr_path(file.path(paths$favor_dir, match_meth, "gds"), favor_pat, sample_chr)
stopifnot(file.exists(sample_gds))
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
  pcs_arg <- list(scores = scores, id_col = pcs_cfg$id_col, cols = pcs_cfg$cols)
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
out_rds <- file.path(paths$pheno_dir, paste0(run_name, "_", null_or(g0("pheno_out_name"), "pheno_covar.rds")))
saveRDS(bundle, out_rds)
write_dataprep_provenance(
  paths$pheno_dir, step = "assemble_pheno_covar",
  script = "00-data-prep/assemble-pheno-covar.R",
  source_alias = basename(pheno_csv), source_path = pheno_csv,
  config_snapshot = list(run = run_name, outcome = outcome_cfg, covariates = g0("covariates"),
                         pcs = pcs_cfg, trait = g0("trait")))
cat(sprintf("Bundle: %d samples, %d covariates [%s] -> %s\n",
            length(bundle$Y), ncol(bundle$X), paste(bundle$covar_names, collapse = ", "), out_rds))
