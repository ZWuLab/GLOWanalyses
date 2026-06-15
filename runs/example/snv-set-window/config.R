# ============================================================================
# runs/example/snv-set-window/config.R - WINDOW-style SNV-set smoke (example)
#
# A region_type = "window" SWAP of the example SNV-set run: same cohort, same
# pre-built B/PI, but the sliding-window scan. Used to SMOKE that the window
# path resolves on the synthetic cohort (one chunk -> a non-empty per-chunk table);
# gene is the fully-validated default. Window geometry is sized for the synthetic
# chromosome spans (chr21 ~0.8 Mb, chr22 ~2.2 Mb).

source("base-config/snv_set_base.R", local = TRUE)
source("runs/example/_cohort.R",      local = TRUE)

region_type <- "window"

# ---- Analysis (a)GDS: the FAVOR-annotated tree ----
gds_dir     <- file.path(data_root, paste0(base_name, "_gds_favor"), "flexible", "gds")
gds_pattern <- favor_gds_pattern

# ---- Pheno + pre-built B/PI (same as the gene example) ----
pheno_path   <- file.path(data_root, "pheno_covar.rds")   # tracked pre-built bundle
b_func       <- NULL
b_model_path <- file.path(data_root, "b_model.rds")
pi_model_dir <- file.path(data_root, "pi_models")

filter_spec <- list(rare_maf_cutoff = 0.5, variant_type = "SNV",
                    min_mac = 1L, min_variants = 1L)
staar_enabled <- FALSE

# ---- Window geometry (sized for the small synthetic spans) ----
window_size    <- 50000L     # 50 kb windows
step_size      <- 25000L     # 25 kb step (50% overlap)
chunk_strategy <- "fixed_mb"
chunk_size_mb  <- 5          # one chunk per chromosome at this scale

primary_test <- "GLOW_Omni"
alpha        <- 0.05
