# ============================================================================
# runs/example/snv-set/config.R - 03-snv-set + 04-summary run for the example
#
# The synthetic-example SNV-set run. Sources base-config/snv_set_base.R then the
# shared cohort block, and points B/PI at the TRACKED pre-built intermediates
# (cohort/b_model.rds + cohort/pi_models/) so 03/04 run STANDALONE on a fresh clone
# (no upstream stages required). region_type = "gene" is the fully-validated
# default; window/coding are documented swaps (see the README) that run through the same scan.
#
# Leaves `symlinked_shared_root` at the base default (NULL = local), so the example
# produces a plain local outputs/ and needs no /project mount.

source("base-config/snv_set_base.R", local = TRUE)
source("runs/example/_cohort.R",      local = TRUE)

# ---- Analysis style + scope ----
region_type <- "gene"          # window / coding are documented region_type swaps
# chroms = c(21L, 22L) comes from the cohort block.

# ---- Analysis (a)GDS: the FAVOR-annotated tree (carries the 16 features) ----
gds_dir     <- file.path(data_root, paste0(base_name, "_gds_favor"), "flexible", "gds")
gds_pattern <- favor_gds_pattern   # chr{chr}_favor.gds

# ---- Pheno bundle: the TRACKED pre-built copy at the cohort root ----
# (build-intermediates.sh copies the 00/assemble-pheno-covar bundle here so 03/04
# run standalone on a fresh clone without the git-ignored 00 lineage dirs.)
pheno_path <- file.path(data_root, "pheno_covar.rds")

# ---- Pre-built B + PI (the TRACKED cohort-root intermediates) ----
b_func       <- NULL                                  # use the trained B model
b_model_path <- file.path(data_root, "b_model.rds")
pi_model_dir <- file.path(data_root, "pi_models")

# ---- Variant filter: GLOW-native full-spectrum (MAF cutoff 0.5) ----
# The synthetic genes hold mostly rare variants + a few common; min_variants = 1
# keeps single-variant sets so every scaffold gene yields a testable region.
filter_spec <- list(rare_maf_cutoff = 0.5, variant_type = "SNV",
                    min_mac = 1L, min_variants = 1L)

# ---- STAAR comparison OFF (keeps the demo STAAR-free + fast) ----
staar_enabled <- FALSE

# ---- Aggregation / significance ----
primary_test <- "GLOW_Omni"
alpha        <- 0.05
top_k        <- 20L            # small cohort -> a short top-hits table
