#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/plink-to-gds.R - convert a cohort's per-chromosome PLINK to GDS
#
# Cohort-agnostic. Wraps GLOWr::plink_to_gds over the PLINK files of `data.plink_dir` and
# writes the GDS files to <output_root>/<dataset>/data/gds/. One call converts one
# chromosome (--chr N), which is one unit of work. The `glow` command starts the units of
# all `shared.chroms`, one after another or as SLURM jobs (`run.mode`). Each unit writes its
# own provenance/config_snapshot_chr<N>.rds.
#
# Usage (from the GLOWanalyses directory):
#   ./glow run data_prep --config <study>.yml --part plink_to_gds     # the usual call: every unit of the part
#   Rscript 00-data-prep/plink-to-gds.R --config <study>.yml --chr 22     # this script, for one unit

suppressMessages(library(GLOWr))
source("00-data-prep/_dataprep_lib.R")
source("lib/glow-stage.R")
st  <- glow_stage_begin("data_prep", "plink_to_gds", options = c(chr = "value"))
pa  <- list(config = st$config_path, opts = st$opts)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

chroms <- if (!is.null(pa$opts$chr)) as.integer(pa$opts$chr) else as.integer(null_or(g0("chroms"), 1:22))
paths  <- stage_cohort_paths(st, g0)
base_name <- paths$base_name
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
  config_snapshot = list(chroms = chroms, plink_pattern = plink_pattern, gds_pattern = gds_pattern),
  unit = if (!is.null(pa$opts$chr)) paste0("chr", pa$opts$chr))
cat("PLINK->GDS complete:", paths$gds_dir, "\n")

# The completion record of this unit.
glow_stage_end(st, outputs = chr_path(paths$gds_dir, gds_pattern, chroms))
