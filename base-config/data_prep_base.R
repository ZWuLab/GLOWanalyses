# ============================================================================
# data_prep_base.R - Documented defaults for the 00-data-prep workflow templates
#
# Cohort-agnostic DEFAULT config for GLOWanalyses/00-data-prep (PLINK->GDS, FAVOR
# annotation, PCs, pheno/covar assembly). A per-run config.R SOURCES this base and
# sets the cohort-specific fields (data_root, base_name, file patterns, the
# pheno/covar mapping) - those NAME real data and so live only in the run config.
#
# DATA-LINEAGE CONVENTION:
#   Data-prep outputs are cohort DATA ASSETS, not run outputs. They live in a
#   user `data_root` as a sibling-suffix tree (<base_name>_gds, _gds_favor, _pcs,
#   _pheno), which the analysis runs then point at via gds_dir/gds_pattern/
#   pheno_path. The templates DERIVE these dirs from data_root + base_name
#   (resolve_cohort_paths in 00-data-prep/_dataprep_lib.R); only file NAME
#   patterns (with a literal {chr} token) live here.

# ---------------------------------------------------------------------------
# Shared-output location — NOT per-run born-shared for 00-data-prep
# ---------------------------------------------------------------------------
# DEFAULT = LOCAL, and the 00-data-prep templates deliberately do NOT consult this
# field: their outputs are cohort DATA ASSETS in the user's `data_root` lineage tree
# (the <base_name>_gds/_gds_favor/_pcs/_pheno siblings, often OUTSIDE the repo), not
# run-org outputs/ leaves, so per-run born-sharing does not apply. Share data-prep
# outputs at the `data_root` level instead (place that root on a shared mount).
# The field is kept here only for base-config parity.
symlinked_shared_root <- NULL

# ---- Cohort lineage anchor (REQUIRED in the run config; no default) ----
# data_root <- ".../processed/<cohort>"   # parent dir holding the sibling tree
# base_name <- "als_landers02_qc_GRCh38"  # the <cohort>_<qc>_<build> identity

# ---- Chromosomes to process ----
chroms <- 1:22

# ---- File-name patterns (literal "{chr}" expands per chromosome) ----
plink_pattern      <- "chr{chr}_hg38"            # PLINK prefix (no extension) in <base>_plink/
gds_pattern        <- "chr{chr}_hg38.gds"        # output in <base>_gds/
favor_gds_pattern  <- "chr{chr}_hg38_favor.gds"  # output in <base>_gds_favor/<match>/gds/
favor_csv_pattern  <- "chr{chr}_hg38_favor.csv"  # output in <base>_gds_favor/<match>/csv/

# ---- plink-to-gds.R (GLOWr::plink_to_gds) ----
plink_to_gds_opts <- list()   # extra args (e.g. list(chr.conv = TRUE))

# ---- annotate-favor.R (GLOWr::annotate_favor) ----
# favor_db        <- ".../FAVOR_annotation/Essential_database_hg38"   # REQUIRED in run config
favor_match_method <- "flexible"   # "exact" (STAAR-compatible) or "flexible" (allele-swap tolerant)
favor_features     <- NULL          # NULL -> annotate_favor() default feature set
favor_use_xsv      <- TRUE
favor_na_handling  <- "keep"

# ---- compute-pcs.R (GLOWr::compute_pcs_gds); reads the aGDS (favor) tree ----
pc_source       <- "favor"    # "favor" (annotated) or "gds" (unannotated) - which tree to read
n_pcs           <- 20L
pc_maf_threshold  <- 0.05
pc_missing_rate   <- 0.05
pc_ld_threshold   <- 0.2
pc_num_thread     <- 2L
pc_seed           <- 42L

# ---- assemble-pheno-covar.R (GLOWr::assemble_pheno_covar) ----
# Reads the sample order from one GDS (sample_gds_chr below), an external pheno
# table, and a PC table; writes <base>_pheno/<run>_pheno_covar.rds.
# REQUIRED in the run config (no cohort-agnostic default):
#   pheno_csv        path to the external phenotype/covariate table
#   pheno_id_col     column in pheno_csv matching the GDS sample.id
#   outcome          list(node= <gds sample.annotation node> OR col= <pheno col>, map= <named vec>)
#   covariates       named list of covariate specs (see assemble_pheno_covar docs)
#   pcs              list(path=, id_col=, cols=)   PC table + its id col + which PCs
trait            <- "binary"
sample_gds_chr   <- 22L        # which chromosome's GDS to read the sample order from
drop_incomplete  <- TRUE
pheno_out_name   <- "pheno_covar.rds"   # prefixed by the run name -> <run>_pheno_covar.rds
