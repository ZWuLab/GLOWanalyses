# ============================================================================
# runs/example/data-prep/config.R - 00-data-prep run for the synthetic example
#
# Drives plink-to-gds -> annotate-favor -> compute-pcs -> assemble-pheno-covar
# over the synthetic example cohort. Sources base-config/data_prep_base.R then the
# shared cohort block (_cohort.R), and sets the data-prep specifics (the outcome /
# covariate / PC mapping for assemble-pheno-covar).
#
# Data-prep outputs are cohort DATA ASSETS in the <data_root> lineage tree (NOT a
# run outputs/ leaf), so this config's only run-org role is to be the --config the
# 00 scripts read; they write into <data_root>/<base_name>_{gds,gds_favor,pcs,pheno}/.

source("base-config/data_prep_base.R", local = TRUE)
source("runs/example/_cohort.R",        local = TRUE)

# ---- PC computation: read the plain (unannotated) GDS tree; few PCs ----
# PCs are genotype-only, so the unannotated example_gds/ tree is the natural source
# (and avoids a SeqArray annotation-attribute merge issue when merging aGDS files).
pc_source       <- "gds"
n_pcs           <- 5L        # small synthetic cohort -> a handful of PCs
pc_maf_threshold <- 0.05      # common-tier variants drive the PCA
pc_missing_rate  <- 0.05
pc_ld_threshold  <- 0.5       # lenient prune (few variants)
pc_seed          <- 42L

# ---- assemble-pheno-covar: outcome + covariates + PCs from the cohort tables ----
sample_gds_chr <- 22L          # read the sample order from chr22's aGDS
trait          <- "binary"

# Outcome from the pheno.csv `outcome` column (already 0/1; identity map).
outcome <- list(col = "outcome", map = c("0" = 0L, "1" = 1L))

# Covariates: sex (already 0/1) + age (numeric) from pheno.csv.
covariates <- list(
  sex = "sex",
  age = "age")

# PCs: the compute-pcs output (a glow_pcs data.frame with a sample.id column and
# PC1..PCn). Match by sample.id; take PC1..PC5 (n_pcs above = 5).
pcs <- list(
  path   = file.path(data_root, paste0(base_name, "_pcs"), "pcs.rds"),
  id_col = "sample.id",
  cols   = paste0("PC", 1:5))

pheno_out_name <- "pheno_covar.rds"   # -> example_pheno/data-prep_pheno_covar.rds
