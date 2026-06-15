#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/compute-pcs.R - principal components from a cohort's GDS tree
#
# Cohort-agnostic. Wraps GLOWr::compute_pcs_gds over the per-chromosome (a)GDS
# files of a cohort lineage tree, writing PCs to <data_root>/<base_name>_pcs/.
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 00-data-prep/compute-pcs.R --config <run>/config.R

suppressMessages(library(GLOWr))
source("00-data-prep/_dataprep_lib.R")

pa  <- parse_dataprep_args()
source(pa$config, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

data_root <- g0("data_root"); base_name <- g0("base_name")
if (is.null(data_root) || is.null(base_name))
  stop("Config must set `data_root` and `base_name`.")
chroms      <- as.integer(null_or(g0("chroms"), 1:22))
paths       <- resolve_cohort_paths(data_root, base_name)

# Which tree to read: annotated (favor) or unannotated (gds).
pc_source   <- match.arg(null_or(g0("pc_source"), "favor"), c("favor", "gds"))
match_meth  <- null_or(g0("favor_match_method"), "flexible")
in_dir      <- if (pc_source == "favor") file.path(paths$favor_dir, match_meth, "gds") else paths$gds_dir
in_pattern  <- if (pc_source == "favor") null_or(g0("favor_gds_pattern"), "chr{chr}_hg38_favor.gds") else null_or(g0("gds_pattern"), "chr{chr}_hg38.gds")
gds_files   <- vapply(chroms, function(c) chr_path(in_dir, in_pattern, c), character(1))
stopifnot(all(file.exists(gds_files)))

dir.create(paths$pcs_dir, recursive = TRUE, showWarnings = FALSE)
out_rds <- file.path(paths$pcs_dir, "pcs.rds")
cat(sprintf("Computing %d PCs over %d GDS (%s tree) -> %s\n",
            as.integer(null_or(g0("n_pcs"), 20L)), length(gds_files), pc_source, out_rds))

pcs <- compute_pcs_gds(
  gds_files     = gds_files,
  n_pcs         = as.integer(null_or(g0("n_pcs"), 20L)),
  maf_threshold = null_or(g0("pc_maf_threshold"), 0.05),
  missing_rate  = null_or(g0("pc_missing_rate"), 0.05),
  ld_threshold  = null_or(g0("pc_ld_threshold"), 0.2),
  output_file   = out_rds,
  num_thread    = as.integer(null_or(g0("pc_num_thread"), 2L)),
  seed          = as.integer(null_or(g0("pc_seed"), 42L)),
  verbose       = 1L)

write.csv(pcs, file.path(paths$pcs_dir, "pcs.csv"), row.names = FALSE)
write_dataprep_provenance(
  paths$pcs_dir, step = "compute_pcs_gds",
  script = "00-data-prep/compute-pcs.R",
  source_alias = paste0(base_name, if (pc_source == "favor") "_gds_favor" else "_gds"),
  source_path = in_dir,
  config_snapshot = list(n_pcs = g0("n_pcs"), maf_threshold = g0("pc_maf_threshold"),
                         missing_rate = g0("pc_missing_rate"), ld_threshold = g0("pc_ld_threshold"),
                         seed = g0("pc_seed"), chroms = chroms, gds_files = gds_files))
cat("PCs written:", out_rds, "(+ pcs.csv, README, provenance/)\n")
