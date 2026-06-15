# ============================================================================
# single_variant_base.R - Documented defaults for the 02-single-variant templates
#
# Cohort-agnostic DEFAULT config for GLOWanalyses/02-single-variant (per-variant
# association scan, in-sample LD scores, genomic-control calibration). A per-run
# config.R SOURCES this base and sets the cohort-specific fields that NAME real
# data (pheno_path, gds_dir, gds_pattern, ld_scores_path) - those live only in
# the run config. Knobs below carry the proven marginal-analysis settings.
#
# RUN-ORG CONVENTION: each 02-single-variant template takes `--config <run>/config.R`,
# DERIVES output_dir = <dir of --config>/outputs and run_name = that dir's name,
# and writes results/plots/logs under that run's outputs/. This is an ANALYSIS
# stage: it READS cohort data assets (the pheno bundle + the cohort GDS), it does
# NOT write into the cohort-data lineage tree (that is 00-data-prep's job).
#
# REQUIRED in the run config (no cohort-agnostic default - they name real data):
#   pheno_path     path to the assemble_pheno_covar bundle (.rds) from 00-data-prep:
#                  a list with $sample_id,$Y,$X,$trait,$covar_names,$n_excluded,
#                  $n_total. fit_null_model is fitted from $X/$Y/$trait/$sample_id.
#   gds_dir        directory holding the per-chromosome cohort (a)GDS files
#   gds_pattern    per-chr GDS file pattern with a literal "{chr}" token
#                  (e.g. "chr{chr}_hg38_favor.gds")
#   ld_scores_path output/input path for the in-sample LD-score table (compute-ld-
#                  scores.R writes it; calibrate.R reads it). LD scores are
#                  GENOTYPE-ONLY, so this is a SHARED location across covariate
#                  runs - point every run's config at the SAME file (do NOT nest
#                  it under a single run's outputs/).

# ---------------------------------------------------------------------------
# Shared-output location — where this run's outputs/ physically live (OPT-IN)
# ---------------------------------------------------------------------------
# DEFAULT = LOCAL (portable; works without a /project shared mount):
#   symlinked_shared_root <- NULL   # outputs/ is a plain local dir under the repo
# OPT IN to born-sharing by overriding in the run's config.R (after it sources this
# base) with a shared-storage root:
#   symlinked_shared_root <- "/project/<lab>/<area>/GLOWpipeline-shared"
#       -> the 02-single-variant templates symlink the run's outputs/ into
#          <root>/<repo-relative-path> (off /home, lab-readable); git stores only
#          the link. On a fresh clone, recreate the symlink before re-running.
symlinked_shared_root <- NULL

# ---- Chromosomes to process (scan + LD; a single-chr run overrides via --chr) ----
chroms <- 1:22

# ---- marginal-scan.R (GLOWr::marginal_scan) ----
use_SPA            <- TRUE     # SPA for a binary trait with rare variants
chunk_size         <- 2000L    # variants per scan chunk
mac_cutoff         <- 1L       # include all non-monomorphic variants (MAC >= 1)
missing_imputation <- "mean"   # genotype missing-value imputation ("mean" or "zero")

# ---- compute-ld-scores.R (GLOWr::compute_ld_scores) ----
ld_window  <- 1e6L   # +/- 1 Mb LD-score window
ld_segment <- 2000L  # core variants per BLAS segment
ld_sample_n <- 0L    # reference-individual subsample size (0 = use all samples)

# ---- calibrate.R (GLOWr::estimate_inflation_factor / calibrate_pvalues) ----
# "ldsc_intercept": LD Score regression intercept (confounding-only; needs LD
# scores; preserves polygenic signal). "lambda_gc": genomic-control lambda.
calibration_method <- "ldsc_intercept"

# NOTE: output_dir is intentionally NOT set here. The templates derive it from the
# --config path (output_dir = <dir of --config>/outputs; run_name = that dir's
# name), so config<->outputs co-location is structural and a run is relocatable.
# Do not add it.
