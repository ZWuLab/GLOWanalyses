# ============================================================================
# runs/example/_cohort.R - the synthetic example cohort definition (shared block)
#
# The SINGLE place that names the synthetic example data. Sourced by every
# per-stage config under runs/example/ (data-prep/, training/, single-variant/,
# snv-set/) AFTER that config sources its base-config, so a stage's config is just
# `source(base); source(_cohort); <a few stage overrides>`. This mirrors the
# research-module run-org pattern (a cohort block + per-stage runs) but for one
# self-contained synthetic cohort that ships in the repo.
#
# All paths are GLOWanalyses-root-relative (run every stage FROM the GLOWanalyses
# directory). The data lives under data-example/cohort/ - the raw inputs + the
# pre-built aGDS + B/PI models are TRACKED (so a clone runs out-of-the-box); the
# per-run GDS/PCs/pheno-bundle/scan/results are generated and git-ignored.
#
# SHARED-STORAGE INDEPENDENCE: every stage's config leaves `symlinked_shared_root`
# at its base default (NULL = local), so the example needs no /project mount.

# ---- Cohort lineage anchor (00-data-prep resolve_cohort_paths) ----
# data_root holds the <base_name>_* sibling tree. base_name = "example" ->
# example_plink/ (raw), example_gds/ + example_gds_favor/ (00 outputs),
# example_pcs/ + example_pheno/ (00 outputs). The pre-built aGDS lives in
# example_gds_favor/flexible/gds/ (built by build-intermediates.sh).
data_root <- "data-example/cohort"
base_name <- "example"

# ---- Chromosomes (the synthetic cohort spans chr21 + chr22) ----
chroms <- c(21L, 22L)

# ---- File-name patterns (literal "{chr}" expands per chromosome) ----
# Override the base-config hg38 defaults for this cohort's simpler names. These
# are the 00-data-prep meanings: plink_pattern -> example_plink/, gds_pattern ->
# example_gds/ (plink-to-gds output), favor_* -> example_gds_favor/ (aGDS).
plink_pattern     <- "chr{chr}"             # example_plink/chr{chr}.{bed,bim,fam}
gds_pattern       <- "chr{chr}.gds"         # example_gds/chr{chr}.gds (plink-to-gds out)
favor_gds_pattern <- "chr{chr}_favor.gds"   # example_gds_favor/flexible/gds/chr{chr}_favor.gds
favor_csv_pattern <- "chr{chr}_favor.csv"   # example_gds_favor/flexible/csv/chr{chr}_favor.csv

# NOTE: the ANALYSIS stages (02-single-variant, 03-snv-set) read the FAVOR aGDS
# tree, so THEIR configs set `gds_dir` + `gds_pattern` to the favor tree (NOT here,
# to avoid clobbering the 00-data-prep gds_pattern above). See those configs.

# ---- Tiny bundled FAVOR database (annotate-favor + PI annotation) ----
favor_db           <- file.path(data_root, "favor-db")
favor_match_method <- "flexible"

# ---- The 16 PI / weighting features carried by the example aGDS ----
# Matches GLOWr:::.default_PI_features() and the trained example PI models.
pi_features <- c(
  "apc_conservation_v2", "apc_epigenetics", "apc_epigenetics_active",
  "apc_epigenetics_repressed", "apc_epigenetics_transcription",
  "apc_protein_function_v3", "apc_local_nucleotide_diversity_v3",
  "apc_mutation_density", "apc_transcription_factor", "apc_mappability",
  "apc_proximity_to_tsstes", "apc_proximity_to_coding_v2", "apc_micro_rna",
  "cadd_phred", "linsight", "fathmm_xf")

# ---- The 3 string coding nodes the coding masks read (GENCODE + MetaSVM) ----
# The coding SNV-set masks resolve plof/missense/... by reading these as string
# annotation sub-nodes (genecode_comprehensive_category /
# genecode_comprehensive_exonic_category / metasvm_pred). The tiny FAVOR DB carries
# them, and annotate_favor now writes them as native string sub-nodes.
coding_nodes <- c("genecode_comprehensive_category",
                  "genecode_comprehensive_exonic_category", "metasvm_pred")

# ---- The full feature set the example aGDS must carry ----
# 00/annotate-favor requests exactly these from the tiny FAVOR DB (16 PI numeric +
# the 3 string coding nodes = the 19 columns the DB provides), so the produced aGDS
# feeds BOTH PI training/prediction and the coding masks. (Left NULL -> the
# annotate_favor() 20-feature default, which this 19-column DB does not fully carry.)
favor_features     <- c(pi_features, coding_nodes)

# ---- Phenotype / covariate table (the cohort's pheno.csv) ----
pheno_csv    <- file.path(data_root, "pheno.csv")
pheno_id_col <- "sample_id"

# ---- Raw known-SNP table (B + PI case training) ----
known_snps_csv <- file.path(data_root, "known-snps.csv")
