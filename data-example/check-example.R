#!/usr/bin/env Rscript
# ============================================================================
# data-example/check-example.R - schema self-test for the synthetic example
#
# Asserts PRESENCE + SCHEMA + NON-EMPTINESS of every stage's artifacts after a
# run-example.sh run. This is the example's CI-able smoke test: it validates the
# template MECHANICS (each stage emitted a schema-valid, non-empty artifact), NOT
# statistical correctness - it deliberately does NOT assert p-value / AUC
# MAGNITUDES (meaningless on synthetic data). Returns a NON-ZERO exit status on
# any missing / wrong-schema / empty artifact.
#
# Usage (from the GLOWanalyses directory, after run-example.sh):
#   conda activate r_env
#   Rscript data-example/check-example.R [output root]     # default: runs

suppressMessages({
  library(SeqArray)
  library(gdsfmt)
  library(data.table)
  library(GLOWr)   # load_B_model / predict_B / load_PI_models
})

# ---- Paths: the derived layout of the example (example.yml) ----
# The dataset's data are under runs/example/data/, the results of the analysis
# `baseline` under runs/example/baseline/. An argument names another output root.
ROOT <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1L] else "runs"
DATA <- file.path(ROOT, "example", "data")
EX   <- file.path(ROOT, "example", "baseline")

# ---- Tiny assertion helper (collects failures; non-zero exit at the end) ----
.fail <- character(0)
ok <- function(cond, msg) {
  status <- isTRUE(cond)
  cat(sprintf("  [%s] %s\n", if (status) "PASS" else "FAIL", msg))
  if (!status) .fail <<- c(.fail, msg)
  invisible(status)
}
file_ok <- function(path, msg = NULL)
  ok(file.exists(path), if (is.null(msg)) sprintf("exists: %s", path) else msg)
# SKIP (not a failure): for artifacts that exist only after a full REBUILD run
# (the GENERATED 00/01 lineage outputs). A default run-example.sh runs only the single-variant
# and SNV-set steps over the pre-built files, so these are absent and intentionally skipped.
# After a rebuild (run-example.sh sets GLOW_EXAMPLE_REBUILT=1 under REBUILD=1) they have to exist.
REBUILT <- identical(Sys.getenv("GLOW_EXAMPLE_REBUILT"), "1")
skip <- function(msg) {
  if (REBUILT) return(invisible(ok(FALSE, paste(msg, "(the run rebuilt the data preparation and the training, so this has to exist)"))))
  cat(sprintf("  [SKIP] %s (generated only; run REBUILD=1 to check)\n", msg))
}

# The 16 PI features the aGDS must carry.
PI_FEATURES <- c(
  "apc_conservation_v2", "apc_epigenetics", "apc_epigenetics_active",
  "apc_epigenetics_repressed", "apc_epigenetics_transcription",
  "apc_protein_function_v3", "apc_local_nucleotide_diversity_v3",
  "apc_mutation_density", "apc_transcription_factor", "apc_mappability",
  "apc_proximity_to_tsstes", "apc_proximity_to_coding_v2", "apc_micro_rna",
  "cadd_phred", "linsight", "fathmm_xf")

# ===========================================================================
cat("=== Stage 00: data-prep ===\n")
# --- Per-chr plain GDS + the FAVOR aGDS carrying the 16 features ---
for (chr in c(21L, 22L)) {
  plain <- file.path(DATA, "gds", sprintf("chr%d.gds", chr))
  agds  <- file.path(DATA, "gds_favor", "flexible", "gds",
                     sprintf("chr%d_favor.gds", chr))
  # Plain GDS is a GENERATED 00 output (not tracked) -> conditional.
  if (file.exists(plain)) ok(TRUE, sprintf("00 plain GDS chr%d", chr))
  else skip(sprintf("00 plain GDS chr%d", chr))
  # aGDS is a TRACKED pre-built asset -> always asserted.
  if (file_ok(agds, sprintf("00 aGDS chr%d", chr))) {
    g <- seqOpen(agds, readonly = TRUE)
    nvar <- length(seqGetData(g, "variant.id"))
    # The aGDS must carry every PI feature as a readable sub-node, non-empty.
    fa_node <- tryCatch(index.gdsn(g, "annotation/info/FunctionalAnnotation"),
                        error = function(e) NULL)
    fa_cols <- if (!is.null(fa_node)) ls.gdsn(fa_node) else character(0)
    have16 <- all(PI_FEATURES %in% fa_cols)
    seqClose(g)
    ok(nvar > 0L, sprintf("00 aGDS chr%d non-empty (%d variants)", chr, nvar))
    ok(have16, sprintf("00 aGDS chr%d carries all 16 PI features", chr))
  }
}
# --- PCs: n_samples x (sample.id + n_pcs) [GENERATED 00 output -> conditional] ---
pcs_path <- file.path(DATA, "pcs", "pcs.rds")
if (!file.exists(pcs_path)) {
  skip("00 PCs (pcs.rds)")
} else if (file_ok(pcs_path, "00 PCs (pcs.rds)")) {
  pcs <- readRDS(pcs_path)
  pc_cols <- grep("^PC", names(pcs), value = TRUE)
  ok(nrow(pcs) > 0L && length(pc_cols) >= 1L && "sample.id" %in% names(pcs),
     sprintf("00 PCs schema: %d samples x %d PCs (+ sample.id)", nrow(pcs), length(pc_cols)))
  ok(!anyNA(as.matrix(pcs[, pc_cols, drop = FALSE])), "00 PCs are finite (no NA)")
}
# --- Pheno bundle: $sample_id / $Y / $X / $trait aligned ---
# Prefer the 00 lineage output; fall back to the tracked cohort-root copy (which is
# what a standalone clone running only 02-04 has).
bundle_path <- file.path(EX, "data_prep", "pheno_covar.rds")
if (file_ok(bundle_path, "00 pheno bundle")) {
  b <- readRDS(bundle_path)
  has_fields <- all(c("sample_id", "Y", "X", "trait") %in% names(b))
  aligned <- has_fields && length(b$sample_id) == length(b$Y) &&
             nrow(b$X) == length(b$Y) && length(b$Y) > 0L
  ok(has_fields, "00 pheno bundle has $sample_id/$Y/$X/$trait")
  ok(aligned, sprintf("00 pheno bundle aligned (%d samples x %d covars)",
                      length(b$Y), if (has_fields) ncol(b$X) else NA))
  ok(b$trait %in% c("binary", "continuous"), sprintf("00 trait valid (%s)", b$trait))
}

# ===========================================================================
cat("=== Stage 01: training (B + PI) ===\n")
# --- B model loads + predicts ---
b_path <- file.path(EX, "training", "b", "b_model.rds")
if (file_ok(b_path, "01 B model (b_model.rds)")) {
  bm <- tryCatch(load_B_model(b_path), error = function(e) NULL)
  ok(!is.null(bm), "01 B model loads")
  if (!is.null(bm)) {
    pred <- tryCatch(predict_B(bm, c(0.001, 0.01, 0.1)), error = function(e) NULL)
    ok(!is.null(pred) && length(pred) == 3L && all(is.finite(pred)),
       "01 B model predicts finite B for 3 target MAFs")
  }
}
# --- PI ensemble: >= 1 model loads ---
pi_dir <- file.path(EX, "training", "pi", "pi_models")
if (ok(dir.exists(pi_dir), "01 PI models dir")) {
  pir <- tryCatch(load_PI_models(pi_dir), error = function(e) NULL)
  ok(!is.null(pir) && pir$n_models >= 1L,
     sprintf("01 PI ensemble loads (%s model(s))",
             if (is.null(pir)) "0" else pir$n_models))
  # the number of models that the training states it trained
  meta <- file.path(pi_dir, "training_metadata.rds")
  if (file_ok(meta, "01 PI training metadata") && !is.null(pir)) {
    stated <- as.integer(readRDS(meta)$n_models)
    ok(identical(as.integer(pir$n_models), stated),
       sprintf("01 PI ensemble holds the %d model(s) that the training states (%d found)", stated, as.integer(pir$n_models)))
  }
}
# --- PI eval table has an AUC column [GENERATED 01 output -> conditional] ---
auc_csv <- file.path(EX, "training", "pi", "pi_eval_auc_per_model.csv")
if (!file.exists(auc_csv)) {
  skip("01 PI eval AUC table")
} else if (file_ok(auc_csv, "01 PI eval AUC table")) {
  auc <- as.data.frame(fread(auc_csv))
  ok(nrow(auc) >= 1L && "auc" %in% names(auc), "01 PI eval table has an `auc` column, non-empty")
}

# ===========================================================================
cat("=== Stage 02: single-variant ===\n")
sv_res <- file.path(EX, "single_variant", "results")
marg <- file.path(sv_res, "single_variant_all.csv")
if (file_ok(marg, "02 single_variant_all.csv")) {
  m <- as.data.frame(fread(marg))
  need_cols <- c("chr", "pos", "Z", "pvalue")
  ok(nrow(m) > 0L && all(need_cols %in% names(m)),
     sprintf("02 single_variant_all schema (chr,pos,...,Z,pvalue); %d variants", nrow(m)))
}
cal <- file.path(sv_res, "single_variant_all_calibrated.csv")
if (file_ok(cal, "02 calibrated CSV")) {
  d <- as.data.frame(fread(cal))
  ok(nrow(d) > 0L && "pvalue" %in% names(d), "02 calibrated CSV non-empty with pvalue")
}

# ===========================================================================
cat("=== Stage 03: SNV-set (gene) ===\n")
shared <- file.path(EX, "snv_set", "gene", "shared")
file_ok(file.path(shared, "01-null_model.rds"), "03 Stage-1 null model")
file_ok(file.path(shared, "01-config_used.rds"), "03 Stage-1 config snapshot")
res_dir <- file.path(EX, "snv_set", "gene", "results")
for (chr in c(21L, 22L)) {
  rc <- file.path(res_dir, sprintf("glow_chr%d.csv", chr))
  if (file_ok(rc, sprintf("03 gene results chr%d", chr))) {
    r <- as.data.frame(fread(rc))
    gene_flat <- all(c("label", "chr", "start", "end", "n_variants") %in% names(r)) &&
                 any(grepl("^GLOW_", names(r)))
    ok(nrow(r) > 0L && gene_flat,
       sprintf("03 chr%d gene-flat schema (label,chr,start,end,...,GLOW_*); %d gene(s)",
               chr, nrow(r)))
  }
}

# ===========================================================================
cat("=== Stage 04: summary ===\n")
agg_dir <- file.path(EX, "snv_set", "gene", "aggregated")
all_csv <- file.path(agg_dir, "glow_results_all.csv")
if (file_ok(all_csv, "04 aggregated table")) {
  a <- as.data.frame(fread(all_csv))
  ok(nrow(a) > 0L && "label" %in% names(a) && any(grepl("^GLOW_", names(a))),
     sprintf("04 aggregated table non-empty with GLOW_* columns (%d genes)", nrow(a)))
}
file_ok(file.path(agg_dir, "glow_top_hits.csv"), "04 top-hits table")
file_ok(file.path(agg_dir, "manhattan_omni.pdf"), "04 Manhattan plot")
ok(length(list.files(agg_dir, pattern = "^qq.*\\.pdf$")) >= 1L, "04 QQ plot(s)")

# ===========================================================================
cat("=== The records of the run ===\n")
for (f in c(file.path(EX, "single_variant", "config-resolved.yml"), file.path(EX, "snv_set", "gene", "config-resolved.yml")))
  file_ok(f, sprintf("resolved configuration: %s", sub(paste0("^", ROOT, "/"), "", f)))
records <- list.files(file.path(EX, "snv_set", "gene", "records"), pattern = "\\.yml$")
ok(all(c("prepare.yml", "scan.chr21.yml", "scan.chr22.yml", "aggregate.yml", "plots.yml") %in% records),
   "completion records of the gene run: prepare, scan chr21, scan chr22, aggregate, plots")

# ===========================================================================
# The state of every unit of the steps that ran, by the rules of the `glow` command: a unit
# is complete when its completion record states the dependencies that example.yml gives now
# and every output of the record exists. This covers every region type of the configuration
# and every chunk of the window scan, and it fails while a job of a unit is still running.
cat("=== The state of every unit ===\n")
source("lib/glow-run.R")
local({
  cfg <- GLOWpipeline::glow_absolute_paths(GLOWpipeline::read_glow_config(
    "example.yml", base = "base-config/base_config.yml", registry = "base-config/settings-registry.yml", root = "."))
  # The units are checked under the output root that was given: the configuration of the example
  # with that root. A folder that is not there, or that this configuration does not own, fails
  # the check. It was skipped before, and a root without its window and coding runs passed.
  cfg$output_root <- normalizePath(ROOT, mustWork = FALSE)
  lay <- GLOWpipeline::glow_layout(cfg)
  if (!dir.exists(EX) || !identical(normalizePath(lay$analysis_dir, mustWork = FALSE), normalizePath(EX))) {
    ok(FALSE, sprintf("the units cannot be checked: %s is not the analysis folder of example.yml under the output root %s", EX, ROOT))
    return(invisible(NULL))
  }
  parts <- GLOWpipeline::glow_read_parts("base-config/parts.yml", "base-config/settings-registry.yml")
  for (step in c(if (REBUILT) c("data_prep", "training"), "single_variant", "snv_set")) {
    region_types <- if (step == "snv_set") cfg$snv_set$region_type else list(NULL)
    for (rt in region_types) for (spec in Filter(function(p) identical(p$step, step), parts)) {
      if (!glow_run_wanted(spec, cfg, rt, lay, parts)$wanted) next
      label <- paste(c(step, rt, spec$part), collapse = " ")
      units <- glow_run_units(spec, cfg, lay, rt, parts)
      if (!length(units)) { ok(FALSE, sprintf("units of %s: none is known (its preparation has not written them)", label)); next }
      states <- vapply(units, function(u) glow_unit_status(cfg, lay, parts, spec, rt, u$unit)$status, character(1))
      open <- vapply(units[!states %in% c("complete", "reused")], function(u) if (is.null(u$unit)) "the part" else u$unit, character(1))
      ok(!length(open), sprintf("units of %s: %d of %d complete%s", label, sum(states %in% c("complete", "reused")), length(units),
                                if (length(open)) paste0(" (not complete: ", paste(utils::head(open, 5L), collapse = ", "), ")") else ""))
    }
  }
})

cat("\n==================================================================\n")
if (length(.fail) == 0L) {
  cat("  ALL CHECKS PASSED.\n")
  cat("==================================================================\n")
  quit(save = "no", status = 0L)
} else {
  cat(sprintf("  %d CHECK(S) FAILED:\n", length(.fail)))
  for (f in .fail) cat("    - ", f, "\n", sep = "")
  cat("==================================================================\n")
  quit(save = "no", status = 1L)
}
