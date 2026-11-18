#!/usr/bin/env Rscript
# ============================================================================
# 02-single-variant/compute-ld-scores.R - in-sample LD scores
#
# Cohort-agnostic step that computes per-variant LD scores from the cohort's own
# genotypes (the SAME (a)GDS the marginal scan uses, so the LD scores cover the
# scan's variants) via GLOWr::compute_ld_scores(), per chromosome, into one
# cached table.
#
# GENOTYPE-ONLY: LD scores are independent of phenotype/covariates. By default
# the table is written into the single_variant folder of the analysis,
# <output_root>/<dataset>/<analysis>/single_variant/ld_scores/ld_scores.rds, and an
# analysis that extends this one reuses it while the settings it depends on are the
# same. `single_variant.ld_scores_path` names another place, for example one file
# that several studies share.
#
# Output: `ld_scores_path` - a data.frame(chr, pos, ref, alt, ld). calibrate.R reads it.
#
# Two modes, as marginal-scan.R has them (one of the two is given):
#   --chr N     ONE chromosome -> ld_chr<N>.rds beside ld_scores_path (the part
#               `ld_scores`); --combine follows.
#   --combine   rbind the ld_chr*.rds beside ld_scores_path in the order of
#               `chroms` -> ld_scores_path (the part `ld_scores_combine`). No GDS
#               is read. Each ld_chr<N>.rds carries the LD settings that made it
#               (ld_window, ld_segment, ld_sample_n, ld_ref_maf_min, gds_pattern);
#               --combine refuses a file whose settings differ from this config's
#               (a stale file of another run), and stops when a requested chromosome
#               has no file.
#
# Reference SNVs (`ld_ref_maf_min`, config; default 0 = every biallelic SNV):
#   above 0, only the SNVs with MAF >= ld_ref_maf_min enter the LD scores, as
#   targets and as neighbors, read from this run's single-variant table
#   (results/single_variant_chr<N>.csv, else single_variant_all.csv), so marginal-scan.R
#   runs first. Restricting the reference SNVs reduces the computation (the
#   all-SNV table is out of reach on whole-genome data) and can change the
#   estimated LDSC intercept: the kept SNVs' LD scores change, and the fit is
#   not scale-free. The threshold travels with the per-chromosome files;
#   validate the resulting calibration on the cohort. A chromosome with no SNV
#   above the threshold stops the script with a message naming it.
#
# Usage (from the GLOWanalyses directory):
#   ./glow run single_variant --config <study>.yml --part ld_scores     # the usual call: every unit of the part
#   Rscript 02-single-variant/compute-ld-scores.R --config <study>.yml --chr 22     # this script, for one unit
#   Rscript 02-single-variant/compute-ld-scores.R --config <study>.yml --combine     # this script, for one unit
# With --combine the script is the part ld_scores_combine. About 10 to 25 minutes genome-wide.

suppressMessages(library(GLOWr))   # GLOWr::compute_ld_scores()
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <path> [--chr N | --combine] ----
source("lib/glow-stage.R")
st <- glow_stage_begin("single_variant", "ld_scores", options = c(chr = "value", combine = "flag"),
                       part_if = c(combine = "ld_scores_combine"))
config_path  <- st$config_path
chr          <- if (is.null(st$opts$chr)) NULL else as.integer(st$opts$chr)
combine_only <- isTRUE(st$opts$combine)
if (!is.null(chr) && combine_only) stop("--chr and --combine are mutually exclusive.")
if (!is.null(chr)) stopifnot(chr >= 1L, chr <= 22L)
# The stem of the result files, "single_variant" (release 0.3.6 wrote "marginal").
sv_stem      <- st$sv_stem
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

# ---- Required cohort fields + LD knobs (settings of the configuration) ----
gds_dir <- g0("gds_dir"); gds_pattern <- g0("gds_pattern"); ld_scores_path <- g0("ld_scores_path")
if (is.null(gds_dir) || is.null(gds_pattern) || is.null(ld_scores_path))
  stop("Config must set `gds_dir`, `gds_pattern`, and `ld_scores_path`.")
ld_window      <- as.integer(null_or(g0("ld_window"), 1e6L))
ld_segment     <- as.integer(null_or(g0("ld_segment"), 2000L))
ld_sample_n    <- as.integer(null_or(g0("ld_sample_n"), 0L))
ld_ref_maf_min <- as.numeric(null_or(g0("ld_ref_maf_min"), 0))
n_cores        <- as.integer(null_or(g0("n_cores"), 1L))   # forked workers over the LD segments (GLOWr 0.3.5)
chroms         <- if (!is.null(chr)) chr else null_or(g0("chroms"), 1:22)
stopifnot(is.finite(ld_ref_maf_min), ld_ref_maf_min >= 0, ld_ref_maf_min < 0.5)
# The settings that determine a chromosome's LD table. --chr N writes them into its
# ld_chr<N>.rds as an attribute and --combine checks them, so a file of another
# configuration in the same directory is refused instead of combined.
ld_settings <- list(ld_window = ld_window, ld_segment = ld_segment, ld_sample_n = ld_sample_n,
                    ld_ref_maf_min = ld_ref_maf_min, gds_pattern = gds_pattern)

gds_path_for <- function(k) file.path(gds_dir, gsub("{chr}", k, gds_pattern, fixed = TRUE))
# A chromosome's own file, beside the shared table: written by --chr N, read by --combine.
per_chr_path <- function(k) file.path(dirname(ld_scores_path), sprintf("ld_chr%s.rds", k))
# This run's single-variant tables (rider: the MAF restriction reads MAF from them).
results_dir  <- file.path(st$output_dir, "results")
# The only output of this step is the LD table at `ld_scores_path` (see the header for its
# default place and for sharing one table among studies).
dir.create(dirname(ld_scores_path), recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# --combine: assemble the per-chromosome files (no GDS read) -> ld_scores_path
# ===========================================================================
if (combine_only) {
  all_chr <- null_or(g0("chroms"), 1:22)
  files   <- per_chr_path(all_chr)
  have    <- file.exists(files)
  if (!all(have))
    stop(sprintf("%d of %d requested chromosomes have no ld_chr<N>.rds beside %s (missing: %s). Every chromosome of `chroms` must have its LD file before --combine.",
                 sum(!have), length(have), ld_scores_path, paste(all_chr[!have], collapse = ", ")))
  # rbind keeps the per-chromosome row names ("21.1", ...), so the table equals the one the
  # all-chromosome call writes, row names included.
  read_checked <- function(f) {
    x <- readRDS(f)
    s <- attr(x, "ld_settings")
    if (is.null(s))
      stop("ld_chr file without LD settings (written by an earlier version?): ", f,
           "\nRecompute it with --chr before combining.")
    if (!identical(s, ld_settings))
      stop("ld_chr file computed with other LD settings than this config: ", f,
           "\n  file  : ", paste(names(s), unlist(s), sep = "=", collapse = ", "),
           "\n  config: ", paste(names(ld_settings), unlist(ld_settings), sep = "=", collapse = ", "),
           "\nRecompute it with --chr, or point ld_scores_path at another directory.")
    attr(x, "ld_settings") <- NULL   # the combined table carries no attribute, so it equals the all-chromosome call's
    x
  }
  ld_all <- do.call(rbind, lapply(files[have], read_checked))   # the order of chroms, as the all-chromosome call
  saveRDS(ld_all, ld_scores_path)
  cat(sprintf("Combined %d chromosome files -> %s : %d variants, mean LD score %.2f\n",
              sum(have), ld_scores_path, nrow(ld_all), mean(ld_all$ld)))
  # the completion record of the part ld_scores_combine
  glow_stage_end(st, outputs = ld_scores_path)
  quit(save = "no", status = 0L)
}

# ---- Reference SNVs for one chromosome: all (default) or MAF >= ld_ref_maf_min ----
# Reads the scan table marginal-scan.R wrote for this run; variant_id restricts both the
# targets and the LD-score neighbors inside compute_ld_scores().
reference_ids <- function(k) {
  if (ld_ref_maf_min <= 0) return(NULL)
  f <- file.path(results_dir, sprintf("%s_chr%s.csv", sv_stem, k))
  if (!file.exists(f)) f <- file.path(results_dir, paste0(sv_stem, "_all.csv"))
  if (!file.exists(f))
    stop("ld_ref_maf_min = ", ld_ref_maf_min, " needs this run's single-variant table ",
         "(results/", sv_stem, "_chr", k, ".csv or ", sv_stem, "_all.csv); run marginal-scan.R first.")
  scan <- as.data.frame(data.table::fread(f, select = c("chr", "variant_id", "MAF"),
                                          colClasses = list(character = "chr")))   # the three columns the selection needs
  # The scan keeps the file's chromosome label ("21" or "chr21"); compare without the prefix
  scan <- scan[sub("^[Cc][Hh][Rr]", "", scan$chr) == as.character(k), , drop = FALSE]
  stopifnot(all(c("variant_id", "MAF") %in% names(scan)))
  ids <- scan$variant_id[is.finite(scan$MAF) & scan$MAF >= ld_ref_maf_min]
  cat(sprintf("chr%-2s  reference SNVs: %d of %d with MAF >= %.3g (from %s)\n",
              k, length(ids), nrow(scan), ld_ref_maf_min, basename(f)))
  as.integer(ids)
}

# ---- Every requested chromosome's GDS must exist before any LD score is computed ----
# A missing file was skipped with a warning, so a wrong path gave a table missing a chromosome
gds_paths <- vapply(chroms, gds_path_for, character(1))
absent    <- gds_paths[!file.exists(gds_paths)]
if (length(absent) > 0L)
  stop("GDS file(s) missing for requested chromosome(s): ", paste(absent, collapse = ", "),
       ". Every chromosome of `chroms` must have its file (gds_dir, gds_pattern).")

# ---- Compute per chromosome (keep cols, incremental saveRDS) ----
res <- list()
t0  <- proc.time()[3]
for (k in chroms) {
  gds_file <- gds_path_for(k)
  tc <- proc.time()[3]
  ids <- reference_ids(k)
  if (!is.null(ids) && length(ids) == 0L)
    stop(sprintf("chr%s: no reference SNVs passed MAF >= %.3g; lower ld_ref_maf_min or drop the chromosome from chroms.",
                 k, ld_ref_maf_min))
  ld <- compute_ld_scores(gds_file, variant_id = ids, window = ld_window,
                          segment = ld_segment, sample_n = ld_sample_n, verbose = 0,
                          n_cores = n_cores)
  res[[as.character(k)]] <- ld[, c("chr", "pos", "ref", "alt", "ld")]
  # --chr N: this chromosome's own file (combine later); default: the shared table, incrementally.
  out_path <- if (!is.null(chr)) per_chr_path(chr) else ld_scores_path
  tbl <- do.call(rbind, res)
  if (!is.null(chr)) attr(tbl, "ld_settings") <- ld_settings   # --chr N: stamp the settings for --combine
  saveRDS(tbl, out_path)   # incremental: a partial run is usable
  cat(sprintf("chr%-2s  n=%6d  mean LD=%.2f  median=%.2f  (%.0fs)\n",
              k, nrow(ld), mean(ld$ld), stats::median(ld$ld),
              proc.time()[3] - tc)); flush.console()
}
ld_all <- do.call(rbind, res)
out_path <- if (!is.null(chr)) per_chr_path(chr) else ld_scores_path
cat(sprintf("\nWrote %s : %d variants, %d chromosomes, mean LD score %.2f, total %.0fs%s\n",
            out_path, nrow(ld_all), length(res), mean(ld_all$ld), proc.time()[3] - t0,
            if (!is.null(chr)) sprintf("\nNext: run --combine once every chromosome has its ld_chr<N>.rds -> %s", ld_scores_path) else ""))

# The completion record of this unit.
glow_stage_end(st, outputs = out_path)
