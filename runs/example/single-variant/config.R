# ============================================================================
# runs/example/single-variant/config.R - 02-single-variant run for the example
#
# Drives marginal-scan -> compute-ld-scores -> calibrate over the synthetic cohort.
# Sources base-config/single_variant_base.R then the shared cohort block. Reads the
# pheno bundle from 00-data-prep and the FAVOR aGDS tree; writes the scan + the
# calibrated CSV into this run's outputs/ (single-variant/outputs/).

source("base-config/single_variant_base.R", local = TRUE)
source("runs/example/_cohort.R",            local = TRUE)

# ---- Analysis (a)GDS: the FAVOR-annotated tree (carries the 16 features) ----
gds_dir     <- file.path(data_root, paste0(base_name, "_gds_favor"), "flexible", "gds")
gds_pattern <- favor_gds_pattern   # chr{chr}_favor.gds

# ---- Pheno bundle: the TRACKED pre-built copy at the cohort root ----
# (build-intermediates.sh copies the 00/assemble-pheno-covar bundle to
# cohort/pheno_covar.rds so this stage runs standalone on a fresh clone. The
# canonical 00 output is example_pheno/data-prep_pheno_covar.rds.)
pheno_path <- file.path(data_root, "pheno_covar.rds")

# ---- LD scores: a SHARED genotype-only asset (kept in the cohort tree) ----
ld_scores_path <- file.path(data_root, paste0(base_name, "_ld"), "ld_scores.rds")

# ---- Scan knobs (small cohort; SPA on for the binary trait) ----
use_SPA            <- TRUE
chunk_size         <- 2000L
mac_cutoff         <- 1L
calibration_method <- "ldsc_intercept"
