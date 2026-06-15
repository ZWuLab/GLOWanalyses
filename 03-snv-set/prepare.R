#!/usr/bin/env Rscript
# ============================================================================
# 03-snv-set/prepare.R - Stage 1: shared prep for any SNV-set GLOW scan
#
# Cohort-agnostic. Runs once per run. Fits the null model, snapshots the
# resolved config, and writes the shared artifacts every Stage-2 task reads.
# The shared Stage-1 setup (dirs, GDS-presence check, pheno load + sample
# alignment + GLOWr null-model fit, PI load + feature resolve, B-source
# dispatch, checksums, base snapshot) is GLOWpipeline::prepare_scan_run(); the
# region_type-specific extra steps run after it and fold into the snapshot:
#   - window : per-chr annotation-median cache + chromosome spans + the genome
#              CHUNK TABLE (GLOWpipeline::define_genome_chunks) + STAAR null.
#   - coding : two STAAR null models (non-SPA fit_null_glm -> STAAR_O_glowG /
#              STAAR_O_native; SPA fit_null_glm_Binary_SPA -> STAAR_B_native)
#              + the essentialdb Annotation_name_catalog + per-chr median cache.
#   - gene   : no Stage-1 extra steps (medians are computed per-chr in run-gene.R).
# GC calibration + use_spa are resolved here for every type.
#
# Outputs (under <output_dir>/shared/):
#   01-null_model.rds          GLOWr null model (GDS sample order)
#   01-config_used.rds         resolved config + git SHA + sessionInfo
#   01-input_checksums.rds     md5 sums for every input file
#   01-annotation_medians.rds  (window/coding) named list: chr -> per-feature medians
#   01-genome_chunks.rds       (window) define_genome_chunks() table
#   01-staar_null_model.rds    (window, if staar_enabled)
#   01-staar_null_model_{spa,nospa}.rds + 01-staar_annotation_catalog.rds (coding)
#
# Usage (from project root):
#   conda activate r_env
#   Rscript 03-snv-set/prepare.R \
#       --config runs/example/config.R

# ---- Setup (dev pattern: GLOWr first, then GLOWpipeline) ----
suppressMessages(library(GLOWr))
suppressMessages(library(GLOWpipeline))
library(SeqArray)

null_or <- function(a, b) if (is.null(a)) b else a

# ---- Args: --config <run>/config.R (required) ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) {
    config_path <- args[i + 1L]; i <- i + 2L
  } else stop("Unknown argument: ", args[i])
}
if (is.na(config_path))
  stop("Missing required --config <run>/config.R")
stopifnot(file.exists(config_path))
cat("Config: ", config_path, "\n", sep = "")

# ---- Source config (VALUES only); sourced at top level = global env ----
source(config_path, local = TRUE)
config <- environment()
chroms <- as.integer(null_or(get0("chroms"), 1:22))
stopifnot(length(chroms) >= 1L, all(chroms >= 1L & chroms <= 22L))

region_type <- match.arg(null_or(get0("region_type"), "gene"),
                         c("gene", "window", "coding"))

# ---- Run identity: output_dir + run_name DERIVED from the --config location ----
if (!is.null(get0("output_dir", inherits = FALSE)))
  message("Note: config set `output_dir`; ignoring it (output_dir is derived ",
          "from the --config path, not declared).")
output_dir <- file.path(dirname(config_path), "outputs")
run_name   <- basename(dirname(config_path))
cat(sprintf("Run: %s [region_type = %s]  ->  output_dir = %s\n",
            run_name, region_type, output_dir))

# ---- Optional config fields -> explicit defaults (required fields fail fast) ----
staar_enabled       <- isTRUE(get0("staar_enabled", ifnotfound = FALSE, inherits = FALSE))
staar_anno_features <- get0("staar_anno_features", ifnotfound = NULL, inherits = FALSE)
use_spa             <- get0("use_spa", ifnotfound = NULL, inherits = FALSE)

# ---- GC calibration: COMPUTE inflation_factor when enabled (else 1.0) ----
calibration      <- get0("calibration", ifnotfound = NULL, inherits = FALSE)
calibration_used <- NULL
if (is.list(calibration) && isTRUE(calibration$enabled)) {
  cal_method <- null_or(calibration$method, "ldsc_intercept")
  if (is.null(calibration$single_variant_scan))
    stop("calibration$enabled = TRUE requires `single_variant_scan` (the scan to ",
         "estimate the inflation factor from; must match this run's covariates).")
  sv_scan <- read.csv(calibration$single_variant_scan, stringsAsFactors = FALSE,
                      colClasses = c(chr = "character"))
  cal_ld  <- if (identical(cal_method, "ldsc_intercept"))
               readRDS(calibration$ld_scores) else NULL
  calibration_used <- estimate_inflation_factor(sv_scan, method = cal_method,
                                                ld_scores = cal_ld)
  inflation_factor <- calibration_used$factor
  if (!is.finite(inflation_factor) || inflation_factor <= 0)
    stop("Computed inflation_factor is not usable: ", inflation_factor)
  message(sprintf("  GC calibration (%s): inflation_factor = %.4f (n = %d) -> z_scale = %.4f (GLOW-only)",
                  cal_method, inflation_factor, calibration_used$n_variants,
                  1 / sqrt(inflation_factor)))
} else {
  inflation_factor <- get0("inflation_factor", ifnotfound = 1.0, inherits = FALSE)
}

# ---- variant_type guard for coding (v1 honours SNV only) ----
if (region_type == "coding") {
  stopifnot(is.list(filter_spec), !is.null(filter_spec$variant_type))
  if (!identical(filter_spec$variant_type, "SNV"))
    stop("Coding masks honour variant_type = 'SNV' only (indels are a deferred ",
         "GLOWr extension); config set '", filter_spec$variant_type, "'.")
}

# ---- Per-type result/sidecar dirs (the shared setup always makes shared/logs/slurm-logs) ----
make_dirs <- switch(region_type,
  gene   = c("results", "logs", "slurm-logs"),
  window = c("results", "logs", "slurm-logs"),
  coding = c("logs", "slurm-logs", "results-glow"))

# ---- Shared Stage-1 setup ----
bundle <- prepare_scan_run(
  config, output_dir = output_dir, chroms = chroms,
  make_dirs = make_dirs, config_path = config_path, verbose = 1)

shared_dir       <- file.path(output_dir, "shared")
null_model       <- bundle$null_model
sample_ids       <- bundle$sample_ids
pi_feats         <- bundle$pi_features
gds_files        <- bundle$gds_files
b_func_set       <- !is.null(bundle$b_func)
b_model_path_set <- !is.null(bundle$base_snapshot$b_model_path)
trait_resolved   <- null_model$trait
bs               <- bundle$base_snapshot

# ---- Common base snapshot fields (every region_type) ----
config_used <- list(
  config_path        = bs$config_path,
  run_name           = run_name,
  region_type        = region_type,
  config_source      = bs$config_source,
  gds_dir            = bs$gds_dir,
  gds_pattern        = bs$gds_pattern,
  pheno_path         = bs$pheno_path,
  chroms             = chroms,
  pi_model_dir       = bs$pi_model_dir,
  b_model_path       = if (b_model_path_set) bs$b_model_path else NULL,
  b_func             = if (b_func_set) bs$b_func else NULL,
  filter_spec        = bs$filter_spec,
  ld_threshold       = bs$ld_threshold,
  mac_threshold      = bs$mac_threshold,
  collapse_method    = bs$collapse_method,
  use_spa            = use_spa,
  inflation_factor   = inflation_factor,
  calibration        = calibration_used,
  pi_features        = bs$pi_features,
  staar_enabled      = staar_enabled,
  staar_anno_features= null_or(staar_anno_features, pi_feats),
  # Output policy (config-authoritative)
  write_flat_table     = TRUE,
  write_evidence       = isTRUE(get0("write_evidence",       ifnotfound = FALSE, inherits = FALSE)),
  write_staar_detail   = isTRUE(get0("write_staar_detail",   ifnotfound = FALSE, inherits = FALSE)),
  write_result_objects = isTRUE(get0("write_result_objects", ifnotfound = FALSE, inherits = FALSE)),
  output_format        = null_or(get0("output_format", ifnotfound = "csv", inherits = FALSE), "csv"),
  # Aggregation / significance
  primary_test           = null_or(get0("primary_test",      ifnotfound = "GLOW_Omni", inherits = FALSE), "GLOW_Omni"),
  primary_test_glow      = null_or(get0("primary_test_glow", ifnotfound = "GLOW_Omni", inherits = FALSE), "GLOW_Omni"),
  alpha                  = null_or(get0("alpha", ifnotfound = 0.05, inherits = FALSE), 0.05),
  staar_genomewide_alpha = null_or(get0("staar_genomewide_alpha", ifnotfound = 2.5e-6, inherits = FALSE), 2.5e-6),
  output_dir         = bs$output_dir,
  trait              = bs$trait,
  glow_test_used_spa = bs$glow_test_used_spa,
  git_sha            = bs$git_sha,
  session_info       = bs$session_info,
  prepared_at        = bs$prepared_at)

input_checksums <- bundle$checksums

# ---- Re-read the aligned pheno (STAAR null fits, when enabled) ----
if (staar_enabled) {
  prep      <- readRDS(config$pheno_path)
  prep_idx  <- match(sample_ids, prep$sample_id)
  Y_aligned <- prep$Y[prep_idx]
  X_aligned <- prep$X[prep_idx, , drop = FALSE]
  if (!requireNamespace("STAAR", quietly = TRUE))
    stop("staar_enabled = TRUE but the STAAR package is not installed.")
}

# ===========================================================================
# region_type-specific extra steps (computed once, folded into the snapshot)
# ===========================================================================
if (region_type == "window") {
  # --- STAAR null model (single, non-SPA omnibus) ---
  if (staar_enabled) {
    cat("Fitting STAAR null model (STAAR::fit_null_glm)...\n")
    staar_df    <- data.frame(Y = Y_aligned, X_aligned, check.names = FALSE)
    family_call <- if (trait_resolved == "binary") binomial(link = "logit") else gaussian()
    staar_null_model <- STAAR::fit_null_glm(Y ~ ., data = staar_df, family = family_call)
    saveRDS(staar_null_model, file.path(shared_dir, "01-staar_null_model.rds"))
    cat("  STAAR null model fitted.\n")
  }

  # --- Per-chr median cache + chromosome spans + the genome chunk table ---
  filter <- do.call(variant_filter, filter_spec)
  cat(sprintf("Computing per-chr medians + spans for chr(s): %s ...\n",
              paste(chroms, collapse = ",")))
  ref_medians_by_chr <- list(); chrom_spans_rows <- list()
  for (chr in chroms) {
    g   <- seqOpen(gds_files[[as.character(chr)]], readonly = TRUE)
    pos <- seqGetData(g, "position")
    chrom_spans_rows[[as.character(chr)]] <- data.frame(
      chr = as.character(chr), start = min(pos), end = max(pos),
      stringsAsFactors = FALSE)
    ref_medians_by_chr[[as.character(chr)]] <- compute_annotation_medians(
      g, pi_feats, filter, sample_id = null_model$sample_id, verbose = 0)
    seqClose(g)
    cat(sprintf("  chr%s: %d variants, span [%d, %d]\n",
                chr, length(pos), min(pos), max(pos)))
  }
  chrom_spans   <- do.call(rbind, chrom_spans_rows)
  genome_chunks <- define_genome_chunks(
    chrom_spans, window_size = window_size, step_size = step_size,
    chunk_strategy = chunk_strategy, chunk_size_mb = chunk_size_mb)
  cat(sprintf("Chunk table: %d chunks; largest est_n_windows = %d\n",
              nrow(genome_chunks), max(genome_chunks$est_n_windows)))

  saveRDS(ref_medians_by_chr, file.path(shared_dir, "01-annotation_medians.rds"))
  saveRDS(genome_chunks,      file.path(shared_dir, "01-genome_chunks.rds"))

  # Window-scan snapshot fields.
  config_used <- c(config_used, list(
    window_size = window_size, step_size = step_size,
    chunk_strategy = chunk_strategy, chunk_size_mb = chunk_size_mb,
    merge_gap   = null_or(get0("merge_gap", ifnotfound = 0L, inherits = FALSE), 0L),
    cmac_cutoff = get0("cmac_cutoff", ifnotfound = NULL, inherits = FALSE),
    staar_null_path = if (staar_enabled) file.path(shared_dir, "01-staar_null_model.rds") else NA_character_,
    staar_version   = if (staar_enabled) as.character(packageVersion("STAAR")) else NA_character_))

} else if (region_type == "coding") {
  # --- categories present + recognized (count + labels for the snapshot/log) ---
  stopifnot(exists("categories", inherits = FALSE), length(categories) >= 1L)
  category_labels <- if (is.list(categories)) names(categories) else categories
  if (is.list(categories) && (is.null(category_labels) || any(!nzchar(category_labels))))
    stop("`categories` given as a list must be NAMED (each entry needs a label).")
  cat(sprintf("Coding categories (%d): %s\n", length(category_labels),
              paste(category_labels, collapse = ", ")))

  staar_modes  <- null_or(get0("staar_modes", ifnotfound = c("spa", "nospa"), inherits = FALSE),
                          c("spa", "nospa"))
  staar_catalog <- NULL
  if (staar_enabled) {
    # Per-mode native-STAAR result dirs.
    for (m in staar_modes)
      dir.create(file.path(output_dir, paste0("results-staar-", m)),
                 recursive = TRUE, showWarnings = FALSE)

    # --- Fit BOTH STAAR null models (id_include + n.pheno = 1L set manually) ---
    stopifnot(trait_resolved == "binary")   # SPA is a binary-trait construct
    staar_df <- data.frame(Y = Y_aligned, X_aligned, check.names = FALSE)
    cat("Fitting STAAR null models (non-SPA fit_null_glm + SPA fit_null_glm_Binary_SPA)...\n")
    staar_null_nospa <- STAAR::fit_null_glm(Y ~ ., data = staar_df,
                                            family = binomial(link = "logit"))
    staar_null_nospa$id_include <- sample_ids; staar_null_nospa$n.pheno <- 1L
    staar_null_spa <- STAAR::fit_null_glm_Binary_SPA(Y ~ ., data = staar_df,
                                                     family = binomial(link = "logit"))
    staar_null_spa$id_include <- sample_ids; staar_null_spa$n.pheno <- 1L
    stopifnot(isTRUE(staar_null_spa$use_SPA), is.null(staar_null_nospa$use_SPA))
    cat(sprintf("  STAAR null models fitted; N = %d\n", length(staar_null_spa$id_include)))

    # --- essentialdb Annotation_name_catalog ---
    staar_catalog <- data.frame(
      name = c("GENCODE.Category", "GENCODE.EXONIC.Category", "GENCODE.Info",
               "MetaSVM", pi_feats),
      dir  = c("/genecode_comprehensive_category",
               "/genecode_comprehensive_exonic_category",
               "/genecode_comprehensive_info", "/metasvm_pred",
               paste0("/", pi_feats)),
      stringsAsFactors = FALSE)

    saveRDS(staar_null_spa,   file.path(shared_dir, "01-staar_null_model_spa.rds"))
    saveRDS(staar_null_nospa, file.path(shared_dir, "01-staar_null_model_nospa.rds"))
    saveRDS(staar_catalog,    file.path(shared_dir, "01-staar_annotation_catalog.rds"))
  }

  # --- Per-chr median cache (+ one-time catalog validation against a live GDS) ---
  Annotation_dir <- null_or(get0("Annotation_dir", ifnotfound = "annotation/info/FunctionalAnnotation",
                                 inherits = FALSE), "annotation/info/FunctionalAnnotation")
  filter <- do.call(variant_filter, filter_spec)
  cat(sprintf("Computing per-chr medians for chr(s): %s ...\n",
              paste(chroms, collapse = ",")))
  ref_medians_by_chr <- list()
  for (chr in chroms) {
    g <- seqOpen(gds_files[[as.character(chr)]], readonly = TRUE)
    if (staar_enabled && identical(chr, chroms[1L])) {
      fa_children <- gdsfmt::ls.gdsn(gdsfmt::index.gdsn(g, Annotation_dir))
      missing_nodes <- setdiff(sub("^/", "", staar_catalog$dir), fa_children)
      if (length(missing_nodes) > 0L) {
        seqClose(g)
        stop("Annotation catalog node(s) missing from chr", chr, " GDS under '",
             Annotation_dir, "': ", paste(missing_nodes, collapse = ", "))
      }
      cat(sprintf("  Catalog validated: all %d nodes resolve in chr%s GDS\n",
                  nrow(staar_catalog), chr))
    }
    ref_medians_by_chr[[as.character(chr)]] <- compute_annotation_medians(
      g, pi_feats, filter, sample_id = null_model$sample_id, verbose = 0)
    pos <- seqGetData(g, "position"); seqClose(g)
    cat(sprintf("  chr%s: %d variants\n", chr, length(pos)))
  }
  saveRDS(ref_medians_by_chr, file.path(shared_dir, "01-annotation_medians.rds"))

  # Coding / native-STAAR snapshot fields.
  staar_omnibus_by_mode <- c(spa = "STAAR-B", nospa = "STAAR-O")
  config_used <- c(config_used, list(
    categories          = categories,
    staar_modes         = staar_modes,
    staar_bakein_glow   = isTRUE(get0("staar_bakein_glow", ifnotfound = TRUE, inherits = FALSE)),
    gene_num_in_array   = null_or(get0("gene_num_in_array", ifnotfound = 500L, inherits = FALSE), 500L),
    Annotation_dir      = Annotation_dir,
    staar_rv_num_cutoff = null_or(get0("staar_rv_num_cutoff", ifnotfound = 2L, inherits = FALSE), 2L),
    geno_missing_imputation = null_or(get0("geno_missing_imputation", ifnotfound = "mean", inherits = FALSE), "mean"),
    staar_write_subtests= isTRUE(get0("staar_write_subtests", ifnotfound = FALSE, inherits = FALSE)),
    staar_omnibus_by_mode = staar_omnibus_by_mode[staar_modes],
    staar_version       = if (staar_enabled) as.character(packageVersion("STAAR")) else NA_character_,
    staar_null_spa_path   = if (staar_enabled) file.path(shared_dir, "01-staar_null_model_spa.rds")   else NA_character_,
    staar_null_nospa_path = if (staar_enabled) file.path(shared_dir, "01-staar_null_model_nospa.rds") else NA_character_,
    staar_catalog_path    = if (staar_enabled) file.path(shared_dir, "01-staar_annotation_catalog.rds") else NA_character_))
}
# (region_type == "gene": no extra steps; medians are computed per-chr in run-gene.R.)

# ---- Persist the snapshot + checksums (null model already written by the shared setup) ----
saveRDS(config_used,     file.path(shared_dir, "01-config_used.rds"))
saveRDS(input_checksums, file.path(shared_dir, "01-input_checksums.rds"))

cat("Stage 1 complete. Shared artifacts at: ", shared_dir, "\n", sep = "")
if (region_type == "window") {
  cat(sprintf("Next: submit the %d-chunk array (slurm/run-window.sh) for run '%s'.\n",
              nrow(genome_chunks), run_name))
} else {
  cat(sprintf("Next: submit the per-chr array (slurm/run-%s.sh) for run '%s'.\n",
              region_type, run_name))
}
