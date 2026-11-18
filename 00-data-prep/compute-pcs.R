#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/compute-pcs.R - principal components from a cohort's GDS tree
#
# Cohort-agnostic. Wraps GLOWr::compute_pcs_gds over the per-chromosome (a)GDS
# files of the dataset, writing the PCs to <output_root>/<dataset>/data/pcs/.
#
# Usage (from the GLOWanalyses directory):
#   ./glow run data_prep --config <study>.yml --part compute_pcs     # the usual call: every unit of the part
#   Rscript 00-data-prep/compute-pcs.R --config <study>.yml     # this script, directly

suppressMessages(library(GLOWr))
source("00-data-prep/_dataprep_lib.R")
source("lib/glow-stage.R")
st  <- glow_stage_begin("data_prep", "compute_pcs")
pa  <- list(config = st$config_path, opts = st$opts)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

chroms      <- as.integer(null_or(g0("chroms"), 1:22))
# Which tree to read: annotated (favor) or unannotated (gds). The folders come from the
# layout (stage_cohort_paths()), never rebuilt here.
pc_source   <- match.arg(null_or(g0("pc_source"), "favor"), c("favor", "gds"))
match_meth  <- null_or(g0("favor_match_method"), "flexible")
paths       <- stage_cohort_paths(st, g0, match_method = match_meth)
base_name   <- paths$base_name
# What is read is the user's folder when the configuration names one (the annotated GDS files by
# shared.analysis_gds_dir with shared.analysis_gds_pattern, the plain ones by
# data_prep.favor_input_gds_dir), and the folder that the run wrote otherwise (paths$read_*).
in_dir      <- if (pc_source == "favor") null_or(paths$read_favor_gds_dir, paths$favor_gds_dir) else null_or(paths$read_gds_dir, paths$gds_dir)
in_pattern  <- if (pc_source == "favor") null_or(paths$read_favor_gds_pattern, null_or(g0("favor_gds_pattern"), "chr{chr}_hg38_favor.gds")) else null_or(g0("gds_pattern"), "chr{chr}_hg38.gds")
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

# The completion record of this unit.
glow_stage_end(st, outputs = c(out_rds, file.path(paths$pcs_dir, "pcs.csv")))
