# ============================================================================
# runs/example/snv-set-coding/config.R - CODING-style SNV-set smoke (example)
#
# A region_type = "coding" SWAP of the example SNV-set run: same cohort, same
# pre-built B/PI, but the per-gene x per-coding-category scan. Used to SMOKE that
# the coding path resolves on the synthetic cohort (the aGDS carries the coding
# nodes genecode_comprehensive_category / _exonic_category / metasvm_pred, which
# the coding masks read). The synthetic gene variants span plof / missense /
# disruptive_missense / synonymous categories. STAAR comparison is OFF (so the
# native STAAR is not a dependency); gene is the fully-validated default.

source("base-config/snv_set_base.R", local = TRUE)
source("runs/example/_cohort.R",      local = TRUE)

region_type <- "coding"

# ---- Analysis (a)GDS: the FAVOR-annotated tree (carries the coding nodes) ----
gds_dir     <- file.path(data_root, paste0(base_name, "_gds_favor"), "flexible", "gds")
gds_pattern <- favor_gds_pattern

# ---- Pheno + pre-built B/PI (same as the gene example) ----
pheno_path   <- file.path(data_root, "pheno_covar.rds")   # tracked pre-built bundle
b_func       <- NULL
b_model_path <- file.path(data_root, "b_model.rds")
pi_model_dir <- file.path(data_root, "pi_models")

# Coding honours variant_type = "SNV" only; min_variants = 1 keeps single-variant
# category sets so a small synthetic gene still yields a testable per-category set.
filter_spec <- list(rare_maf_cutoff = 0.5, variant_type = "SNV",
                    min_mac = 1L, min_variants = 1L)

# ---- Coding categories present in the synthetic gene variants ----
# (The generator seeds stopgain/stoploss/splicing -> plof; nonsynonymous(D/T) ->
# missense / disruptive_missense; synonymous SNV -> synonymous.)
categories <- c("plof", "missense", "disruptive_missense", "synonymous")

staar_enabled     <- FALSE   # GLOW-only coding (no native-STAAR comparator)
staar_bakein_glow <- FALSE   # do not also emit STAAR_O_glowG (needs a STAAR null)

primary_test      <- "GLOW_Omni"
primary_test_glow <- "GLOW_Omni"
alpha             <- 0.05
