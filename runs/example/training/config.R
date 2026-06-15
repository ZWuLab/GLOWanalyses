# ============================================================================
# runs/example/training/config.R - 01-training run for the synthetic example
#
# Drives B (prepare -> estimate) and PI (prepare -> annotate -> train -> evaluate)
# on the synthetic cohort's known-snps.csv + tiny FAVOR DB + the pre-built aGDS
# control tree. Sources base-config/training_base.R then the shared cohort block.
#
# The B model + PI ensemble land in this run's outputs/ (training/outputs/). The
# build-intermediates.sh step ALSO copies them up to the cohort root
# (cohort/b_model.rds, cohort/pi_models/) as the TRACKED pre-built intermediates an
# SNV-set run points at, so 03/04 can run standalone on a fresh clone.

source("base-config/training_base.R", local = TRUE)
source("runs/example/_cohort.R",       local = TRUE)

# ---- B + PI raw known-SNP table (the cohort's known-snps.csv) ----
b_raw_path  <- known_snps_csv
pi_raw_path <- known_snps_csv

# ---- FAVOR features to annotate: the 16 numeric PI features for the PI-case CSV ----
# PI-case annotation + the eval feature set need only the numeric PI features (the
# coding nodes the cohort block adds for the analysis aGDS are not used here), so
# override the cohort block's full 19-column set with exactly the 16.
favor_features <- pi_features

# ---- B training ----
b_trait_type   <- "binary"
b_training_trait <- "binary"
b_method       <- "auto"          # the known-SNP table carries MAF/BETA/P/N
# Relax the default QC min_sample_size (the synthetic N = 3000 passes, but keep
# explicit so a smaller synthetic N would still train).
b_qc_filters   <- list(min_sample_size = 100, max_pvalue = 1, min_maf = 0,
                       remove_duplicates = TRUE, duplicate_priority = "N",
                       remove_na = TRUE, remove_na_cols = c("MAF", "P", "N"))

# ---- PI training ----
# Tiny ensemble for the demo (enough to exercise load_PI_models + evaluate).
pi_n_models           <- 5L
pi_model_type         <- "LASSO"
pi_controls_multiplier <- 2       # 2x controls per case (cases are few)
pi_max_controls       <- 5000L    # the whole cohort has < this many variants
pi_eval_max_controls  <- 2000L
pi_chromosomes        <- c(21L, 22L)
pi_random_seed        <- 42L
pi_models_subdir      <- "pi_models"
# pi_features inherits the cohort block's 16-feature list.

# ---- PI control source: the pre-built FAVOR-annotated aGDS tree ----
# Built by build-intermediates.sh (the aGDS builder). The control loader reads the
# FunctionalAnnotation node from each chr's aGDS (carries the 16 PI features).
pi_control_source <- file.path(data_root, paste0(base_name, "_gds_favor"),
                               "flexible", "gds")
pi_control_format <- "gds"
