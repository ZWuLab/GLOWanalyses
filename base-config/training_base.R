# ============================================================================
# training_base.R - Documented defaults for the 01-training workflow templates
#
# Cohort-agnostic DEFAULT config for GLOWanalyses/01-training (B effect-size-
# distribution estimation + PI prior-inclusion-probability ensemble training).
# A per-run config.R SOURCES this base and sets the cohort-specific fields (the
# raw known-SNP table, the FAVOR database, the cohort control tree, and the
# column mappings / exclusions for that cohort's data) - those NAME real data and
# so live only in the run config.
#
# RUN-ORGANIZATION CONVENTION (mirrors snv_set_base.R):
#   - A RUN is a directory runs/<name>/ holding a config.R and a sibling outputs/.
#   - config.R SOURCES this base, then sets the cohort paths + any overrides.
#   - The --config path passed to every step fixes WHERE outputs live:
#     output_dir = <dir of --config>/outputs ; run_name = that dir's name.
#   - Do NOT set `output_dir` here or in a run config; it is DERIVED.
#
# STEP CHAINING (each step reads the previous step's output in the same run):
#   B  : prepare-b-training-data.R -> b_training_data.rds
#        estimate-b.R              -> b_model.rds (+ b_values.csv)
#   PI : prepare-pi-case.R   -> pi_case_data.csv
#        annotate-pi-case.R  -> pi_case_data_annotated.csv
#        train-pi.R          -> <pi_models_subdir>/model_*.rds
#        evaluate-pi.R       -> pi_eval_*.csv + pi_eval_roc.pdf
#   CONTROL annotations are produced UPSTREAM (00-data-prep/annotate-favor.R or
#   system.file("scripts/annotate_favor_hpc.sh", package = "GLOWr")); a run points train/eval
#   at that already-annotated tree via `pi_control_source`.

# ===========================================================================
# Shared-output location — where this run's outputs/ physically live (OPT-IN)
# ===========================================================================
# DEFAULT = LOCAL (portable; works without a /project shared mount):
#   symlinked_shared_root <- NULL   # outputs/ is a plain local dir under the repo
# OPT IN to born-sharing by overriding in the run's config.R (after it sources this
# base) with a shared-storage root:
#   symlinked_shared_root <- "/project/<lab>/<area>/GLOWpipeline-shared"
#       -> the 01-training templates symlink the run's outputs/ into
#          <root>/<repo-relative-path> (off /home, lab-readable); git stores only
#          the link. On a fresh clone, recreate the symlink before re-running.
symlinked_shared_root <- NULL

# ===========================================================================
# B (effect-size distribution) - prepare-b-training-data.R + estimate-b.R
# ===========================================================================

# ---- Raw known-SNP table (REQUIRED in the run config; no default) ----
# b_raw_path <- ".../ALS-known-SNPs-raw.xlsx"  # known trait-associated SNPs
#   The GLOWr-bundled example is system.file("extdata/ALS-known-SNPs-raw.xlsx", package = "GLOWr").

# ---- prepare-b-training-data.R (GLOWr::prepare_B_training_data) ----
b_trait_type       <- "binary"   # "binary" | "continuous" | "mixed"
b_format           <- "auto"     # file-format hint for the loader
# b_column_mapping : NAMED list mapping the standard names to the raw columns,
#   e.g. list(BETA="Beta_numeric", MAF="RAF", P="P-value", N="Sample Size", ...).
#   Data-specific -> set in the run config. NULL -> loader auto-detect.
b_column_mapping   <- NULL
# b_filter_include / b_filter_exclude : NAMED lists keyed by a raw column, each
#   value a vector of values/patterns to keep / drop (e.g. drop non-standard
#   reported traits or a specific first author). Data-specific -> run config.
b_filter_include   <- NULL
b_filter_exclude   <- NULL
b_filter_pattern   <- TRUE       # treat filter values as substring patterns (vs exact)
b_filter_autosomes <- TRUE       # keep autosomes (chr 1-22) only
# b_qc_filters : NULL -> prepare_B_training_data() built-in QC defaults
#   (min_sample_size=500, max_pvalue=1, min_maf=0, remove_duplicates=TRUE,
#    duplicate_priority="N", remove_na=TRUE, remove_na_cols=c("MAF","P","N")).
#   Set a list to override; a stricter run sets e.g. min_sample_size=1000.
b_qc_filters       <- NULL

# ---- estimate-b.R (GLOWr::train_B_model / save_B_model / predict_B) ----
# b_training_rds : NULL -> <output_dir>/b_training_data.rds (this run's prepare step)
b_training_rds     <- NULL
b_training_trait   <- "binary"   # "binary" | "continuous" | "mixed"
b_method           <- "auto"     # "beta" | "pvalue" | "both" | "auto"
# Which columns IN the training table feed the fit. Defaults match the
# prepare_B_training_data() output schema. Set a *_col to NULL to drop that
# input (e.g. b_p_col=NULL, b_n_col=NULL -> beta-only fit).
b_maf_col          <- "MAF"
b_beta_col         <- "BETA"
b_p_col            <- "P"
b_n_col            <- "N"
# b_target_maf : OPTIONAL. NULL -> skip the B(MAF) table. Either a numeric vector
#   of target MAFs, OR a list(from=, to=, n=) grid spec -> b_values.csv.
b_target_maf       <- NULL
b_target_trait     <- NULL       # NULL -> use b_training_trait
b_target_case_prop <- NULL       # binary-only; NULL for continuous targets

# ===========================================================================
# PI (prior inclusion probability) - prepare/annotate/train/evaluate
# ===========================================================================

# ---- Raw known-SNP table (REQUIRED in the run config; no default) ----
# pi_raw_path <- ".../ALS-known-SNPs-raw.xlsx"   # same family of raw inputs as B

# ---- prepare-pi-case.R (GLOWr::prepare_PI_case_data) ----
pi_format             <- "auto"  # file-format hint for the loader
# pi_column_mapping : NAMED list mapping standard names to raw columns,
#   e.g. list(N="Sample Size", BETA="Beta_numeric"). NULL -> auto-detect.
pi_column_mapping     <- NULL
pi_filter_include     <- NULL    # NAMED list of keep filters (data-specific)
pi_filter_exclude     <- NULL    # NAMED list of drop filters (data-specific)
pi_exclude_authors    <- NULL    # character vector of first-authors to drop (e.g. "Nicolas A")
pi_exclude_chromosomes <- NULL   # character vector of chromosomes to drop (e.g. "X")
# pi_qc_filters : NULL -> prepare_PI_case_data() built-in QC defaults
#   (remove_duplicates=TRUE, duplicate_priority="N", filter_autosomes=TRUE,
#    filter_indels=FALSE, remove_na=TRUE, remove_na_cols=c("rsID","CHR","POS")).
pi_qc_filters         <- NULL

# ---- annotate-pi-case.R (GLOWr::annotate_favor on the case csv) ----
# favor_db : REQUIRED in the run config (names a real FAVOR DB on disk):
#   favor_db <- ".../FAVOR_annotation/Essential_database_hg38"
# pi_case_csv : NULL -> <output_dir>/pi_case_data.csv (this run's prepare step)
pi_case_csv           <- NULL
favor_match_method    <- "flexible"  # "exact" (STAAR-compatible) | "flexible" (allele-swap tolerant)
favor_features        <- NULL        # NULL -> annotate_favor() default feature set
favor_use_xsv         <- TRUE        # fast xsv-based pre-filter
favor_na_handling     <- "keep"      # "keep" | "remove" | "impute"

# ---- train-pi.R (GLOWr::train_PI_models) ----
# pi_case_annotated_csv : NULL -> <output_dir>/pi_case_data_annotated.csv
pi_case_annotated_csv <- NULL
# pi_control_source : REQUIRED in the run config. The cohort's FAVOR-annotated
#   CONTROL tree (a directory of per-chromosome aGDS/CSV), built UPSTREAM by
#   00-data-prep/annotate-favor.R or system.file("scripts/annotate_favor_hpc.sh", package = "GLOWr").
#   pi_control_source <- ".../<cohort>_gds_favor/<match>/gds"
pi_n_models           <- 100L     # ensemble size
pi_controls_multiplier <- 1       # controls_per_model = n_cases * multiplier
pi_max_controls       <- 50000L   # cap on controls loaded (NULL = all); memory knob
pi_model_type         <- "LASSO"  # "GLM" | "LASSO" (LASSO matches legacy, smaller files)
pi_chromosomes        <- 1:22     # chromosomes to sample controls from
pi_features           <- NULL     # NULL -> GLOWr:::.default_PI_features()
pi_random_seed        <- 42L
pi_models_subdir      <- "pi_models"  # models land in <output_dir>/<this> (the SNV-set pi_model_dir)

# ---- evaluate-pi.R (GLOWr::evaluate_PI_models / plot_PI_roc) ----
# pi_models_dir : NULL -> <output_dir>/<pi_models_subdir> (this run's trained ensemble)
pi_models_dir         <- NULL
pi_control_format     <- "auto"   # "auto" | "gds" | "csv" - control file format for eval
pi_eval_max_controls  <- 5000L    # controls loaded for evaluation (NULL = all); lower = faster/noisier
pi_n_eval_controls    <- NULL     # subsample controls for the ROC (NULL = use pi_eval_max_controls)
