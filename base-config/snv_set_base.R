# ============================================================================
# snv_set_base.R - Documented defaults for the SNV-set GLOW workflow templates
#
# This is the cohort-agnostic, documented DEFAULT config for the GLOWanalyses
# 03-snv-set / 04-summary templates. It carries every parameter the unified
# GLOW scan reads, with the safe/recommended default and an inline
# comment for each. It sets NO cohort paths (gds_dir, pheno_path, pi_model_dir,
# b_model_path/b_func) - those name real data and so belong only in a per-run
# config (runs/<name>/config.R), which sources this base and overrides what the
# run varies.
#
# RUN-ORGANIZATION CONVENTION:
#   - A RUN is a directory runs/<name>/ holding a config.R and its sibling
#     outputs/.
#   - config.R SOURCES this base, then sets the cohort paths + any overrides.
#   - The --config path passed to every stage fixes WHERE the outputs live:
#     output_dir = <dir of --config>/outputs ; run_name = that dir's name.
#   - Do NOT set `output_dir` here or in a run config; it is DERIVED, so a run
#     is relocatable and config<->outputs co-location is structural.
#
# region_type SELECTS the analysis style; the thin per-type entry script
# (03-snv-set/run-{gene,window,coding}.R) sets it and parses the unit selector
# (--chr for gene/coding, --chunk for window). The window/coding-only fields
# below are ignored when region_type is something else.

# ---------------------------------------------------------------------------
# region_type - SELECTS the analysis style (gene | window | coding)
# ---------------------------------------------------------------------------
# - "gene"   : test each gene's rare-variant set (one filter, many genes).
# - "window" : sliding-window scan (one filter, windows tiled per chunk).
# - "coding" : per-gene x per-coding-category sets (one filter PER category).
# The per-type entry script sets this explicitly, but a run config may also pin
# it so a bare run-* call matches the run's intent.
region_type <- "gene"

# ---------------------------------------------------------------------------
# Chromosomes to scan
# ---------------------------------------------------------------------------
# Full genome-wide run uses 1:22. A run that overrides this (e.g. chroms <- 22L)
# is a single-chromosome smoke: Stage 1 builds only these chromosomes' caches /
# chunk table, so the array shrinks accordingly.
chroms <- 1:22

# ---------------------------------------------------------------------------
# Cohort paths - NOT SET HERE (a run config must supply them)
# ---------------------------------------------------------------------------
# Required per-run (named for real data):
#   gds_dir      directory of per-chromosome (a)GDS files
#   gds_pattern  filename pattern with the literal token {chr}, e.g.
#                "chr{chr}.gds" (the scripts substitute the chromosome number)
#   pheno_path   .rds with a prepared pheno/covar object: a list with
#                $sample_id (character), $Y (response), $X (covariate matrix),
#                and $trait ("binary" | "continuous") - see GLOWr::fit_null_model
#   pi_model_dir directory of trained PI ensemble models (GLOWr::get_PI*)
# B source (set EXACTLY ONE per run):
#   b_func       a function(maf) -> per-variant B (closed form), OR
#   b_model_path .rds of a trained B model (GLOWr::get_B*)

# ---------------------------------------------------------------------------
# Variant filter
# ---------------------------------------------------------------------------
# filter_spec is passed unchanged to GLOWr::variant_filter (do.call). Fields:
#   rare_maf_cutoff  MAF axis. 0.5 = GLOW-native full-spectrum weighting +
#                    collapsing; a sequencing-style run uses 0.01.
#   variant_type     "SNV" (the only v1-supported value for coding; indels are
#                    a deferred GLOWr extension).
#   min_mac          drop variants monomorphic in the cohort (cohort MAC = 0);
#                    1L matches the GLOWr default and is kept explicit.
#   min_variants     1L keeps single-variant sets (glow_test reduces them to the
#                    single-variant SNV test; the is_single flag marks them).
filter_spec <- list(rare_maf_cutoff = 0.5, variant_type = "SNV",
                    min_mac = 1L, min_variants = 1L)

# ---------------------------------------------------------------------------
# LD pruning / collapsing thresholds (scan knobs)
# ---------------------------------------------------------------------------
ld_threshold    <- 0.95   # prune one of any pair with r^2 above this
mac_threshold   <- 11L    # collapse ultra-rare variants below this cohort MAC
collapse_method <- "mean" # how a collapsed unit's genotype is summarised

# ---------------------------------------------------------------------------
# PI / annotation weighting features
# ---------------------------------------------------------------------------
# pi_features is the set of annotation columns GLOW's PI (and, if STAAR is on,
# STAAR's annotation weights) read from the (a)GDS. NULL falls back to
# GLOWr:::.default_PI_features(). A run typically sets the concrete list that
# matches its trained PI models and the annotations present in its GDS.
pi_features <- NULL

# ---------------------------------------------------------------------------
# SPA option for GLOW's single-variant score statistics (default = auto)
# ---------------------------------------------------------------------------
# Selects the per-unit Z that every GLOW family test combines:
#   NULL  = auto: SPA for binary traits, standard for continuous (the default).
#   TRUE  = force SPA (binary only; a non-binary trait errors).
#   FALSE = force standard score stats for any trait.
# SPA is recommended for rare variants / unbalanced case-control; standard is
# acceptable (and cheaper) for common-variant, balanced, large-n analyses.
# GLOW-only: STAAR runs its own (always non-SPA) path and is unaffected.
use_spa <- NULL

# ---------------------------------------------------------------------------
# Shared-output location — where this run's outputs/ physically live
# ---------------------------------------------------------------------------
# outputs/ ALWAYS appears in the repo; this field only decides real-dir vs SYMLINK,
# and is read by GLOWpipeline::prepare_scan_run (which every SNV-set template routes
# through). For the USER-FACING GLOWanalyses layer the DEFAULT is LOCAL (portable —
# works without a /project shared mount), set explicitly just below:
#
#   symlinked_shared_root <- NULL   # (the default here)
#       outputs/ is a plain local directory under the repo (no symlink).
#
# OPT IN to born-sharing for a particular run by overriding the field in that run's
# config.R (AFTER it sources this base) with a shared-storage root:
#
#   symlinked_shared_root <- "/project/<lab>/<area>/GLOWpipeline-shared"
#       outputs/ becomes a symlink into <root>/<repo-relative-path>; the bytes live
#       there (off /home, lab-readable) and git stores only the link. On a fresh
#       clone, recreate the symlink before re-running.
#
# The resolved choice is recorded in the run snapshot (01-config_used.rds).
symlinked_shared_root <- NULL

# ---------------------------------------------------------------------------
# Genomic-control (GC) calibration of GLOW (default OFF)
# ---------------------------------------------------------------------------
# When calibration$enabled, Stage 1 COMPUTES a single-variant inflation_factor
# via GLOWr::estimate_inflation_factor() from the declared sources (NOT a
# hardcoded value; the snapshot records the provenance), and Stage 2 scales
# every GLOW marginal Z by z_scale = 1/sqrt(inflation_factor) BEFORE GLOW
# combines them. GLOW-ONLY: it never touches any STAAR column (STAAR consumes
# the genotype G, not GLOW's Z). Absent/disabled -> inflation_factor = 1.0 (an
# exact no-op). A directly-set `inflation_factor` is a manual fallback.
#   method "ldsc_intercept" isolates confounding from polygenicity (needs
#   ld_scores); "lambda_gc" needs no ld_scores. single_variant_scan MUST match
#   this run's covariates.
calibration <- list(
  enabled             = FALSE,
  method              = "ldsc_intercept",  # or "lambda_gc"
  single_variant_scan = NULL,              # path to a single-variant scan CSV
  ld_scores           = NULL               # path to an LD-scores .rds
)

# ---------------------------------------------------------------------------
# STAAR comparison (optional; default ON for window/coding, off for gene)
# ---------------------------------------------------------------------------
# staar_enabled toggles the STAAR comparison: when TRUE, Stage 1 fits the STAAR null
# model(s) and Stage 2 builds a staar_context so each set also carries a STAAR
# comparator. STAAR is a hard dependency only when this is TRUE.
#   staar_anno_features  STAAR's annotation-weight columns; NULL falls back to
#                        pi_features (the fair, same-features setting).
staar_enabled       <- FALSE
staar_anno_features <- NULL

# ---------------------------------------------------------------------------
# WINDOW-only fields (ignored unless region_type == "window")
# ---------------------------------------------------------------------------
# window_size / step_size define the sliding geometry; chunk_strategy /
# chunk_size_mb control how Stage 1 splits the genome into SLURM-array chunks.
window_size    <- 100000L  # 100 kb window
step_size      <- 10000L   # 10 kb step
chunk_strategy <- "fixed_mb"
chunk_size_mb  <- 5        # STAARpipeline-proven default
# Native-STAARpipeline sliding-window comparison (optional; run-staar-window.R):
# when TRUE (and staar_enabled), prepare.R also fits both native STAAR null
# models + the annotation catalog, so STAARpipeline's own engine can be run on
# exactly the GLOW window grid, once per mode in staar_modes (fields below).
staar_native   <- FALSE
# Annotation weights in the NATIVE engine. The native path consumes annotation
# values raw (no imputation; only missing CADD is zeroed), so on an aGDS with
# incomplete annotation coverage an NA becomes an NaN weight and crashes the
# SKAT eigendecomposition. Set FALSE there: the native tests then use the pure
# beta(MAF) weights and native STAAR-O equals the CCT of the same six base
# tests as the embedded ACAT_O. (GLOW's embedded comparator is unaffected
# either way: it median-imputes annotations before STAAR sees them.)
staar_native_use_annotation <- TRUE

# ---------------------------------------------------------------------------
# CODING-only fields (ignored unless region_type == "coding")
# ---------------------------------------------------------------------------
# categories is either a character vector of builtin coding-category names OR a
# NAMED list mapping each label to a builtin name / custom DNF annotation
# clauses (build_scan_regions resolves both). The 7 standard STAARpipeline
# coding categories:
categories <- c("plof", "plof_ds", "missense", "disruptive_missense",
                "synonymous", "ptv", "ptv_ds")
# Native-STAARpipeline comparison fields (the optional, supported comparison
# stage; default-driven). Read when staar_enabled: by the coding
# Gene_Centric_Coding stage (run-staar-coding.R), and by the window native
# stage (run-staar-window.R) when staar_native = TRUE above.
staar_modes             <- c("spa", "nospa")  # native STAAR runs once per mode
staar_bakein_glow       <- TRUE   # GLOW also emits STAAR_O_glowG on its G
gene_num_in_array       <- 500L   # genes per native-STAAR job (data-size tuned)
Annotation_dir          <- "annotation/info/FunctionalAnnotation"  # GDS node path
staar_rv_num_cutoff     <- 2L
geno_missing_imputation <- "mean"
staar_write_subtests    <- FALSE  # native STAAR: also save the full STAAR grid

# ---------------------------------------------------------------------------
# Output policy (config-authoritative; the flat table is always written)
# ---------------------------------------------------------------------------
# The flat result table is the load-bearing output; the three sidecars are
# default-off inspection artifacts (written per unit only when the flag is on).
write_flat_table     <- TRUE   # always (glow_output_policy enforces TRUE)
write_evidence       <- FALSE  # per-unit evidence packets (SNV-level traceability)
write_staar_detail   <- FALSE  # per-unit STAAR detail grid (needs a staar_context)
write_result_objects <- FALSE  # per-unit raw glow_test_result S3 objects
output_format        <- "csv"  # "csv" | "fst"

# ---------------------------------------------------------------------------
# Aggregation / significance (Stage 3, 04-summary)
# ---------------------------------------------------------------------------
# primary_test drives locus merging, ranking, significance, and the Manhattan/QQ
# headline. gene/window read `primary_test`; coding reads `primary_test_glow`
# (both default to GLOW_Omni) - set the one matching this run's region_type.
primary_test           <- "GLOW_Omni"  # gene / window
primary_test_glow      <- "GLOW_Omni"  # coding
alpha                  <- 0.05
merge_gap              <- 0L      # window: bp gap; 0 = only overlapping windows merge
cmac_cutoff            <- NULL    # window: optional STAARpipeline-style cMAC filter
top_k                  <- 100L    # gene: top-K hits to write
staar_genomewide_alpha <- 2.5e-6  # coding: STAARpipeline genome-wide convention

# NOTE: output_dir is intentionally NOT set here. The scripts derive it from the
# --config path (see the header). Do not add it.
