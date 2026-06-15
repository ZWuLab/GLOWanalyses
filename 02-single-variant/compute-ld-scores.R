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
# Usage (from project root; ~10-25 min genome-wide):
#   conda activate r_env
#   Rscript 02-single-variant/compute-ld-scores.R \
#       --config <run>/config.R [--chr 22]

suppressMessages(library(GLOWr))   # GLOWr::compute_ld_scores()
null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <path> [--chr N] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; chr <- NULL
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--chr" && i < length(args)) { chr <- as.integer(args[i + 1L]); i <- i + 2L }
  else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))
if (!is.null(chr)) stopifnot(chr >= 1L, chr <= 22L)

# ---- Config sourcing (the proven idiom) ----
source(config_path, local = TRUE)
cfg <- environment()
g0  <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

# ---- Required cohort fields + LD knobs (knobs from base) ----
gds_dir <- g0("gds_dir"); gds_pattern <- g0("gds_pattern"); ld_scores_path <- g0("ld_scores_path")
if (is.null(gds_dir) || is.null(gds_pattern) || is.null(ld_scores_path))
  stop("Config must set `gds_dir`, `gds_pattern`, and `ld_scores_path`.")
ld_window   <- as.integer(null_or(g0("ld_window"), 1e6L))
ld_segment  <- as.integer(null_or(g0("ld_segment"), 2000L))
ld_sample_n <- as.integer(null_or(g0("ld_sample_n"), 0L))
chroms      <- if (!is.null(chr)) chr else null_or(g0("chroms"), 1:22)

gds_path_for <- function(k) file.path(gds_dir, gsub("{chr}", k, gds_pattern, fixed = TRUE))
# Shared-output note: this step does NOT consult `symlinked_shared_root` and does NOT
# write into the run's outputs/. Its only output is the GENOTYPE-ONLY LD table at
# `ld_scores_path`, a SHARED asset across covariate runs (see header). Per-run
# born-sharing therefore does not apply here; share `ld_scores_path` itself (or its
# parent) by placing it on a shared mount if it must be lab-visible. Cf. 00-data-prep,
# which is likewise data_root-level.
dir.create(dirname(ld_scores_path), recursive = TRUE, showWarnings = FALSE)

# ---- Compute per chromosome (keep cols, incremental saveRDS) ----
res <- list()
t0  <- proc.time()[3]
for (k in chroms) {
  gds_file <- gds_path_for(k)
  if (!file.exists(gds_file)) {
    cat(sprintf("chr%-2s  WARNING: %s not found, skipping\n", k, gds_file)); next
  }
  tc <- proc.time()[3]
  ld <- compute_ld_scores(gds_file, window = ld_window, segment = ld_segment,
                          sample_n = ld_sample_n, verbose = 0)
  res[[as.character(k)]] <- ld[, c("chr", "pos", "ref", "alt", "ld")]
  saveRDS(do.call(rbind, res), ld_scores_path)   # incremental: a partial run is usable
  cat(sprintf("chr%-2s  n=%6d  mean LD=%.2f  median=%.2f  (%.0fs)\n",
              k, nrow(ld), mean(ld$ld), stats::median(ld$ld),
              proc.time()[3] - tc)); flush.console()
}
if (length(res) == 0L) stop("No GDS files found for the requested chromosomes; nothing computed.")
ld_all <- do.call(rbind, res)
cat(sprintf("\nWrote %s : %d variants, %d chromosomes, mean LD score %.2f, total %.0fs\n",
            ld_scores_path, nrow(ld_all), length(res), mean(ld_all$ld), proc.time()[3] - t0))
