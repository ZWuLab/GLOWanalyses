#!/usr/bin/env Rscript
# ============================================================================
# 02-single-variant/compute-ld-scores.R - in-sample LD scores
#
# Cohort-agnostic step that computes per-variant LD scores from the cohort's own
# genotypes (the SAME (a)GDS the marginal scan uses, so the LD scores cover the
# scan's variants) via GLOWr::compute_ld_scores(), per chromosome, into one
# cached table.
#
# GENOTYPE-ONLY: LD scores are independent of phenotype/covariates, so this is
# computed ONCE and the one table SERVES EVERY covariate run. Accordingly the
# output `ld_scores_path` is a SHARED location, NOT a per-run outputs/ path -
# point every covariate run's config at the same file. (This script still uses
# the run-org --config entry for a uniform interface; it does NOT write into the
# run's outputs/, because the LD table is a shared genotype asset.)
#
# Output: `ld_scores_path` (config field) - a data.frame(chr, pos, ref, alt, ld),
#         saved incrementally (a partial run is usable). calibrate.R reads it.
#
# Three modes, as marginal-scan.R has them:
#   (default)   loop `chroms` -> ld_scores_path, saved after every chromosome.
#   --chr N     ONE chromosome -> ld_chr<N>.rds beside ld_scores_path (an array
#               task or a cloud task owns one chromosome); --combine follows.
#   --combine   rbind the ld_chr*.rds beside ld_scores_path in the order of
#               `chroms` (as the all-chromosome call) -> ld_scores_path. No GDS
#               is read. Each ld_chr<N>.rds carries the LD settings that made it
#               (ld_window, ld_segment, ld_sample_n, ld_ref_maf_min, gds_pattern);
#               --combine refuses a file whose settings differ from this config's
#               (a stale file of another run), and skips a chromosome without a
#               file with a warning.
#
# Reference SNVs (`ld_ref_maf_min`, config; default 0 = every biallelic SNV):
#   above 0, only the SNVs with MAF >= ld_ref_maf_min enter the LD scores, as
#   targets and as neighbors, read from this run's single-variant table
#   (results/marginal_chr<N>.csv, else marginal_all.csv), so marginal-scan.R
#   runs first. Restricting the reference SNVs reduces the computation (the
#   all-SNV table is out of reach on whole-genome data) and can change the
#   estimated LDSC intercept: the kept SNVs' LD scores change, and the fit is
#   not scale-free. The threshold travels with the per-chromosome files;
#   validate the resulting calibration on the cohort. A chromosome with no SNV
#   above the threshold stops the script with a message naming it.
#
# Usage (from project root; ~10-25 min genome-wide):
#   conda activate r_env
#   Rscript 02-single-variant/compute-ld-scores.R \
#       --config <run>/config.R [--chr 22 | --combine]

suppressMessages(library(GLOWr))   # GLOWr::compute_ld_scores()
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <path> [--chr N | --combine] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; chr <- NULL; combine_only <- FALSE
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--chr" && i < length(args)) { chr <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--combine") { combine_only <- TRUE; i <- i + 1L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))
if (!is.null(chr) && combine_only) stop("--chr and --combine are mutually exclusive.")
if (!is.null(chr)) stopifnot(chr >= 1L, chr <= 22L)

# ---- Config sourcing (the proven idiom) ----
source(config_path, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

# ---- Required cohort fields + LD knobs (knobs from base) ----
gds_dir <- g0("gds_dir"); gds_pattern <- g0("gds_pattern"); ld_scores_path <- g0("ld_scores_path")
if (is.null(gds_dir) || is.null(gds_pattern) || is.null(ld_scores_path))
  stop("Config must set `gds_dir`, `gds_pattern`, and `ld_scores_path`.")
ld_window      <- as.integer(null_or(g0("ld_window"), 1e6L))
ld_segment     <- as.integer(null_or(g0("ld_segment"), 2000L))
ld_sample_n    <- as.integer(null_or(g0("ld_sample_n"), 0L))
ld_ref_maf_min <- as.numeric(null_or(g0("ld_ref_maf_min"), 0))
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
results_dir  <- file.path(dirname(config_path), "outputs", "results")
# Shared-output note: this step does NOT consult `symlinked_shared_root` and does NOT
# write into the run's outputs/. Its only output is the GENOTYPE-ONLY LD table at
# `ld_scores_path`, a SHARED asset across covariate runs (see header). Per-run
# born-sharing therefore does not apply here; share `ld_scores_path` itself (or its
# parent) by placing it on a shared mount if it must be lab-visible. Cf. 00-data-prep,
# which is likewise data_root-level.
dir.create(dirname(ld_scores_path), recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# --combine: assemble the per-chromosome files (no GDS read) -> ld_scores_path
# ===========================================================================
if (combine_only) {
  all_chr <- null_or(g0("chroms"), 1:22)
  files   <- per_chr_path(all_chr)
  have    <- file.exists(files)
  if (!any(have)) stop("No per-chromosome LD files (ld_chr*.rds) beside ", ld_scores_path, " to combine.")
  if (!all(have)) cat(sprintf("WARNING: %d of %d chromosomes have no ld_chr<N>.rds (missing: %s); combining the rest\n",
                              sum(!have), length(have), paste(all_chr[!have], collapse = ", ")))
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
  quit(save = "no", status = 0L)
}

# ---- Reference SNVs for one chromosome: all (default) or MAF >= ld_ref_maf_min ----
# Reads the scan table marginal-scan.R wrote for this run; variant_id restricts both the
# targets and the LD-score neighbors inside compute_ld_scores().
reference_ids <- function(k) {
  if (ld_ref_maf_min <= 0) return(NULL)
  f <- file.path(results_dir, sprintf("marginal_chr%s.csv", k))
  if (!file.exists(f)) f <- file.path(results_dir, "marginal_all.csv")
  if (!file.exists(f))
    stop("ld_ref_maf_min = ", ld_ref_maf_min, " needs this run's single-variant table ",
         "(results/marginal_chr", k, ".csv or marginal_all.csv); run marginal-scan.R first.")
  scan <- read.csv(f, colClasses = list(chr = "character"), check.names = FALSE)
  scan <- scan[scan$chr == as.character(k), , drop = FALSE]
  stopifnot(all(c("variant_id", "MAF") %in% names(scan)))
  ids <- scan$variant_id[is.finite(scan$MAF) & scan$MAF >= ld_ref_maf_min]
  cat(sprintf("chr%-2s  reference SNVs: %d of %d with MAF >= %.3g (from %s)\n",
              k, length(ids), nrow(scan), ld_ref_maf_min, basename(f)))
  as.integer(ids)
}

# ---- Compute per chromosome (keep cols, incremental saveRDS) ----
res <- list()
t0  <- proc.time()[3]
for (k in chroms) {
  gds_file <- gds_path_for(k)
  if (!file.exists(gds_file)) {
    cat(sprintf("chr%-2s  WARNING: %s not found, skipping\n", k, gds_file)); next
  }
  tc <- proc.time()[3]
  ids <- reference_ids(k)
  if (!is.null(ids) && length(ids) == 0L)
    stop(sprintf("chr%s: no reference SNVs passed MAF >= %.3g; lower ld_ref_maf_min or drop the chromosome from chroms.",
                 k, ld_ref_maf_min))
  ld <- compute_ld_scores(gds_file, variant_id = ids, window = ld_window,
                          segment = ld_segment, sample_n = ld_sample_n, verbose = 0)
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
if (length(res) == 0L) stop("No GDS files found for the requested chromosomes; nothing computed.")
ld_all <- do.call(rbind, res)
out_path <- if (!is.null(chr)) per_chr_path(chr) else ld_scores_path
cat(sprintf("\nWrote %s : %d variants, %d chromosomes, mean LD score %.2f, total %.0fs%s\n",
            out_path, nrow(ld_all), length(res), mean(ld_all$ld), proc.time()[3] - t0,
            if (!is.null(chr)) sprintf("\nNext: run --combine once every chromosome has its ld_chr<N>.rds -> %s", ld_scores_path) else ""))
