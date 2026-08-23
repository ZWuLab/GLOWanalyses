#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/annotate-favor.R - annotate a cohort's GDS with FAVOR -> aGDS
#
# Cohort-agnostic. Wraps GLOWr::annotate_favor over <base_name>_gds/ ->
# <base_name>_gds_favor/<match_method>/{gds,csv}/. Runs all `chroms` by default,
# or one chromosome via `--chr N`. FAVOR annotation is the heaviest 00 step, so at
# whole-genome scale run it as a per-chromosome SLURM array via
# slurm/run-annotate-favor.sh (which calls this script with --chr). Each --chr
# task writes its own provenance/config_snapshot_chr<N>.rds.
#
# Usage (from the GLOWanalyses directory):
#   conda activate r_env
#   Rscript 00-data-prep/annotate-favor.R --config <run>/config.R [--chr 22]

suppressMessages(library(GLOWr))
source("00-data-prep/_dataprep_lib.R")

pa  <- parse_dataprep_args()
source(pa$config, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

data_root <- g0("data_root"); base_name <- g0("base_name")
if (is.null(data_root) || is.null(base_name)) stop("Config must set `data_root` and `base_name`.")
favor_db <- g0("favor_db")
if (is.null(favor_db) || !dir.exists(favor_db)) stop("Config `favor_db` must be an existing FAVOR DB dir.")
chroms     <- if (!is.null(pa$opts$chr)) as.integer(pa$opts$chr) else as.integer(null_or(g0("chroms"), 1:22))
match_meth <- null_or(g0("favor_match_method"), "flexible")
features   <- g0("favor_features")   # NULL -> annotate_favor() default set
paths      <- resolve_cohort_paths(data_root, base_name, match_method = match_meth)
gds_pattern   <- null_or(g0("gds_pattern"),        "chr{chr}_hg38.gds")
fav_gds_pat   <- null_or(g0("favor_gds_pattern"),  "chr{chr}_hg38_favor.gds")
fav_csv_pat   <- null_or(g0("favor_csv_pattern"),  "chr{chr}_hg38_favor.csv")

dir.create(paths$favor_gds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$favor_csv_dir, recursive = TRUE, showWarnings = FALSE)
for (chr in chroms) {
  in_gds  <- chr_path(paths$gds_dir,       gds_pattern, chr)
  out_gds <- chr_path(paths$favor_gds_dir, fav_gds_pat, chr)
  out_csv <- chr_path(paths$favor_csv_dir, fav_csv_pat, chr)
  stopifnot(file.exists(in_gds))
  cat(sprintf("chr%d [%s]: %s -> %s\n", chr, match_meth, in_gds, out_gds))
  annot_args <- list(
    variants       = in_gds,
    favor_db_path  = favor_db,
    output_csv     = out_csv,
    output_agds    = out_gds,
    match_method   = match_meth,
    use_xsv        = isTRUE(null_or(g0("favor_use_xsv"), TRUE)),
    na_handling    = null_or(g0("favor_na_handling"), "keep"),
    verbose        = 1)
  # favor_features = NULL means "annotate_favor()'s default feature set": omit the
  # argument so the package default applies (the default is the full FAVOR
  # Essential DB annotation content -- scannable by every built-in category).
  if (!is.null(features)) annot_args$features <- features
  do.call(annotate_favor, annot_args)
}

write_dataprep_provenance(
  paths$favor_dir, step = "annotate_favor",
  script = "00-data-prep/annotate-favor.R",
  source_alias = paste0(base_name, "_gds"), source_path = paths$gds_dir,
  notes = sprintf("match_method = %s; FAVOR DB = %s", match_meth, favor_db),
  config_snapshot = list(chroms = chroms, match_method = match_meth, favor_db = favor_db,
                         features = features),
  unit = if (!is.null(pa$opts$chr)) paste0("chr", pa$opts$chr))
cat("FAVOR annotation complete:", paths$favor_dir, "\n")
