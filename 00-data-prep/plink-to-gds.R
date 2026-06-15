#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/plink-to-gds.R - convert a cohort's per-chromosome PLINK to GDS
#
# Cohort-agnostic. Wraps GLOWr::plink_to_gds over <base_name>_plink/ -> <base_name>_gds/.
# Runs all `chroms` by default, or one chromosome via `--chr N` (SLURM array).
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 00-data-prep/plink-to-gds.R --config <run>/config.R [--chr 22]

suppressMessages(library(GLOWr))
source("00-data-prep/_dataprep_lib.R")

pa  <- parse_dataprep_args()
source(pa$config, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

data_root <- g0("data_root"); base_name <- g0("base_name")
if (is.null(data_root) || is.null(base_name)) stop("Config must set `data_root` and `base_name`.")
chroms <- if (!is.null(pa$opts$chr)) as.integer(pa$opts$chr) else as.integer(null_or(g0("chroms"), 1:22))
paths  <- resolve_cohort_paths(data_root, base_name)
plink_pattern <- null_or(g0("plink_pattern"), "chr{chr}_hg38")
gds_pattern   <- null_or(g0("gds_pattern"),   "chr{chr}_hg38.gds")
extra <- null_or(g0("plink_to_gds_opts"), list())

dir.create(paths$gds_dir, recursive = TRUE, showWarnings = FALSE)
for (chr in chroms) {
  prefix <- chr_path(paths$plink_dir, plink_pattern, chr)
  out    <- chr_path(paths$gds_dir,   gds_pattern,   chr)
  cat(sprintf("chr%d: %s -> %s\n", chr, prefix, out))
  do.call(plink_to_gds, c(list(plink_prefix = prefix, output_gds = out, verbose = 1), extra))
}

write_dataprep_provenance(
  paths$gds_dir, step = "plink_to_gds",
  script = "00-data-prep/plink-to-gds.R",
  source_alias = paste0(base_name, "_plink"), source_path = paths$plink_dir,
  config_snapshot = list(chroms = chroms, plink_pattern = plink_pattern, gds_pattern = gds_pattern))
cat("PLINK->GDS complete:", paths$gds_dir, "\n")
