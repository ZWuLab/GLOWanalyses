#!/usr/bin/env Rscript
# ============================================================================
# 00-data-prep/annotate-favor.R - annotate a cohort's GDS with FAVOR -> aGDS
#
# Cohort-agnostic. Wraps GLOWr::annotate_favor over the GDS files of the dataset
# (data/gds/, or the folder of `data_prep.favor_input_gds_dir`) and writes the annotated
# files to data/gds_favor/<match_method>/{gds,csv}/. One call annotates one chromosome
# (--chr N), which is one unit of work. FAVOR annotation is the heaviest step of the data
# preparation: the `glow` command starts the units of all `shared.chroms`, one after another
# or as SLURM jobs (`run.mode`). Each unit writes its own
# provenance/config_snapshot_chr<N>.rds.
#
# Usage (from the GLOWanalyses directory):
#   ./glow run data_prep --config <study>.yml --part annotate_favor     # the usual call: every unit of the part
#   Rscript 00-data-prep/annotate-favor.R --config <study>.yml --chr 22     # this script, for one unit

suppressMessages(library(GLOWr))
source("00-data-prep/_dataprep_lib.R")
source("lib/glow-stage.R")
st  <- glow_stage_begin("data_prep", "annotate_favor", options = c(chr = "value"))
pa  <- list(config = st$config_path, opts = st$opts)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

favor_db <- g0("favor_db")
if (is.null(favor_db) || !dir.exists(favor_db)) stop("Config `favor_db` must be an existing FAVOR DB dir.")
chroms     <- as.character(if (!is.null(pa$opts$chr)) pa$opts$chr else null_or(g0("chroms"), 1:22))
match_meth <- null_or(g0("favor_match_method"), "flexible")
features   <- g0("favor_features")   # NULL -> annotate_favor() default set
db_format  <- null_or(g0("favor_db_format"), "auto")   # "auto" | "csv" | "parquet"
fav_release <- g0("favor_release")  # NULL -> recorded as "unknown" in the aGDS
rsid_pol   <- null_or(g0("favor_rsid_policy"), "require")  # "require" | "record"
# The input GDS tree and the annotated tree, from the layout: the GDS files of the dataset
# (data/gds/, or the folder of `data_prep.favor_input_gds_dir`) and data/gds_favor/<match>/{gds,csv}/.
paths      <- stage_cohort_paths(st, g0, match_method = match_meth)
base_name  <- paths$base_name
gds_pattern   <- null_or(g0("gds_pattern"),        "chr{chr}_hg38.gds")
fav_gds_pat   <- null_or(g0("favor_gds_pattern"),  "chr{chr}_hg38_favor.gds")
fav_csv_pat   <- null_or(g0("favor_csv_pattern"),  "chr{chr}_hg38_favor.csv")

dir.create(paths$favor_gds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$favor_csv_dir, recursive = TRUE, showWarnings = FALSE)
for (chr in chroms) {
  in_gds  <- chr_path(paths$favor_input_dir, gds_pattern, chr)
  out_gds <- chr_path(paths$favor_gds_dir, fav_gds_pat, chr)
  out_csv <- chr_path(paths$favor_csv_dir, fav_csv_pat, chr)
  stopifnot(file.exists(in_gds))
  cat(sprintf("chr%s [%s]: %s -> %s\n", chr, match_meth, in_gds, out_gds))
  annot_args <- list(
    variants       = in_gds,
    favor_db_path  = favor_db,
    output_csv     = out_csv,
    output_agds    = out_gds,
    match_method   = match_meth,
    use_xsv        = isTRUE(null_or(g0("favor_use_xsv"), TRUE)),
    na_handling    = null_or(g0("favor_na_handling"), "keep"),
    favor_db_format = db_format,
    favor_release  = fav_release,
    rsid_policy    = rsid_pol,
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
  source_alias = paste0(base_name, "_gds"), source_path = paths$favor_input_dir,
  notes = sprintf("match_method = %s; rsid_policy = %s; FAVOR DB = %s (%s)", match_meth, rsid_pol, favor_db, db_format),
  config_snapshot = list(chroms = chroms, match_method = match_meth, rsid_policy = rsid_pol,
                         favor_db = favor_db, favor_db_format = db_format, favor_release = fav_release,
                         features = features,
                         # the resolved paths, so the record says what was read and written
                         favor_input_dir = paths$favor_input_dir,
                         favor_gds_dir = paths$favor_gds_dir, favor_csv_dir = paths$favor_csv_dir),
  unit = if (!is.null(pa$opts$chr)) paste0("chr", pa$opts$chr))
cat("FAVOR annotation complete:", paths$favor_dir, "\n")

# The completion record of this unit.
# Both files are written for every chromosome, so both are outputs of the unit.
glow_stage_end(st, outputs = c(chr_path(paths$favor_gds_dir, fav_gds_pat, chroms),
                               chr_path(paths$favor_csv_dir, fav_csv_pat, chroms)))
