#!/usr/bin/env Rscript
# ============================================================================
# 03-snv-set/run-staar-window.R - Native STAARpipeline sliding-window engine,
# per genome chunk, unit-matched to the GLOW window grid
#
# Runs STAARpipeline's OWN sliding-window analysis (Sliding_Window(type =
# "single"), i.e. the fully native extraction + STAAR/STAAR_Binary_SPA path)
# on EXACTLY the windows the GLOW window scan tests: per chunk it regenerates
# the grid with the same two lines build_scan_regions() uses
# (GLOWr::define_regions_window(chr, chunk_start, span_end, window_size,
# step_size), then keep start <= chunk_end), so every output row joins the
# GLOW table 1:1 by `label`. STAARpipeline's own window generator is NOT used
# because it hard-codes a half-window step (Sliding_Window_Multiple), which
# cannot reproduce a 100 kb / 10 kb grid.
#
# Run ONCE PER MODE (--spa {on,off}); the mode selects the Stage-1 native null
# model (Sliding_Window auto-detects use_SPA) and the output dir:
#   --spa off -> 01-staar_null_model_nospa.rds -> results-staar-nospa/
#                (STAAR-O + ACAT-O + the SKAT/Burden/ACAT-V family grid)
#   --spa on  -> 01-staar_null_model_spa.rds   -> results-staar-spa/
#                (STAAR-B + the Burden family grid)
# SPA_p_filter = FALSE always (the marginal-p pre-filter degenerates small
# sets; same Phase-0 lesson as the coding engine).
#
# The per-family omnibus columns are extracted POSITIONALLY from the native
# return (Sliding_Window names only Chr/Start/End/#SNV/cMAC and the final
# omnibus columns; each family block of length B ends with that family's CCT).
# NOTE ON THE COLUMN NAMES: what is taken is that family's CCT over its
# annotation-weighted variants, i.e. STAAR-S/B/A, not the bare base test. With
# annotation weighting OFF (the `staar_native_use_annotation = FALSE` case) the
# block holds one un-annotated p-value plus its CCT, so the CCT equals the base
# test and the SKAT_*/Burden_*/ACATV_* names are exact. With annotation
# weighting ON they are the STAAR-S/B/A omnibi and the names understate them;
# read them accordingly, or take each block's FIRST element for the base test.
# The row width is asserted against the expected layout (nospa: 5 + 6B + 2;
# spa: 5 + 2B + 1) and the named columns are cross-checked, so a STAAR version
# change breaks loudly, never silently.
#
# Output (under <output_dir>/): results-staar-<mode>/staar_window_<chunk>.csv
# (curated flat table; a 0-row schema frame when the chunk yields nothing) and
# logs/log_staar_window_<mode>_<chunk>.txt.
#
# Usage (from the GLOWanalyses directory; prepare.R must have run with
# staar_native = TRUE so the native null models + catalog exist):
#   conda activate r_env
#   Rscript 03-snv-set/run-staar-window.R \
#       --config runs/example/config.R \
#       --chunk 1 --spa off [--max-windows 5]

suppressMessages({
  library(gdsfmt); library(SeqArray)
  library(STAAR); library(STAARpipeline)
  library(GLOWr)   # define_regions_window: the shared window-grid generator
})

null_or <- function(a, b) if (is.null(a)) b else a
t_start <- Sys.time()

# ---- Args: --config <path> --chunk <id> --spa {on,off} [--max-windows K] ----
args <- commandArgs(trailingOnly = TRUE)
config_path <- NA_character_; chunk_id <- NULL; spa_on <- NA; max_windows <- NA_integer_
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--chunk" && i < length(args)) { chunk_id <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--max-windows" && i < length(args)) { max_windows <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--spa" && i < length(args)) {
    v <- tolower(args[i + 1L])
    spa_on <- if (v %in% c("on", "true", "1")) TRUE
      else if (v %in% c("off", "false", "0")) FALSE
      else stop("--spa must be on|off (got '", args[i + 1L], "')")
    i <- i + 2L
  } else stop("Unknown argument: ", args[i])
}
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
if (is.na(spa_on)) stop("Missing required --spa {on,off}")
stopifnot(file.exists(config_path), !is.null(chunk_id), chunk_id >= 1L)
spa_mode <- if (spa_on) "spa" else "nospa"

# ---- Run identity + Stage-1 artifacts ----
output_dir <- file.path(dirname(config_path), "outputs")
run_name   <- basename(dirname(config_path))
shared_dir <- file.path(output_dir, "shared")
results_staar_dir <- file.path(output_dir, paste0("results-staar-", spa_mode))
logs_dir   <- file.path(output_dir, "logs")
for (d in c(results_staar_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

null_path <- file.path(shared_dir, if (spa_on) "01-staar_null_model_spa.rds"
                                   else        "01-staar_null_model_nospa.rds")
need <- c(file.path(shared_dir, c("01-config_used.rds", "01-genome_chunks.rds",
                                  "01-staar_annotation_catalog.rds")), null_path)
if (any(!file.exists(need)))
  stop("Missing Stage-1 artifacts (run prepare.R with staar_native = TRUE). Missing:\n  ",
       paste(need[!file.exists(need)], collapse = "\n  "))
config        <- readRDS(file.path(shared_dir, "01-config_used.rds"))
genome_chunks <- readRDS(file.path(shared_dir, "01-genome_chunks.rds"))
obj_nullmodel <- readRDS(null_path)          # use_SPA auto-detected by Sliding_Window
staar_catalog <- readRDS(file.path(shared_dir, "01-staar_annotation_catalog.rds"))
# Annotation weighting is optional for the NATIVE path: unlike GLOW's embedded
# comparator (which median-imputes missing annotations), the native engine
# consumes annotation values raw, and an NA becomes an NaN weight that crashes
# STAAR's SKAT eigendecomposition ("eig_sym(): decomposition failed"). On an
# aGDS with incomplete annotation coverage set staar_native_use_annotation to
# FALSE: the native tests then use the pure beta(MAF) weights, and the native
# STAAR-O equals the CCT of the same six base tests as the embedded ACAT_O.
use_anno <- isTRUE(null_or(config$staar_native_use_annotation, TRUE))
Annotation_name <- if (use_anno) null_or(config$staar_anno_features, config$pi_features) else NULL

# ---- Logging ----
log_con <- file(file.path(logs_dir, sprintf("log_staar_window_%s_%04d.txt", spa_mode, chunk_id)),
                open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE, after = FALSE)
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = ""); writeLines(line, log_con); flush(log_con)
}

# ---- Resolve THIS chunk; regenerate ITS windows (the GLOW grid, exactly) ----
chunk <- genome_chunks[genome_chunks$chunk_id == chunk_id, , drop = FALSE]
if (nrow(chunk) != 1L)
  stop(sprintf("chunk_id %d not found in the chunk table (%d chunks).",
               chunk_id, nrow(genome_chunks)))
chr <- as.character(chunk$chr)
# Same two lines as GLOWpipeline::build_scan_regions() (region_type "window"):
# generate to span_end so boundary-spanning windows are kept, then keep only the
# windows whose start this chunk owns.
regions <- GLOWr::define_regions_window(chr, chunk$chunk_start, chunk$span_end,
                                        config$window_size, config$step_size)
regions <- regions[regions$start <= chunk$chunk_end, , drop = FALSE]
rownames(regions) <- NULL
if (!is.na(max_windows) && max_windows >= 1L && nrow(regions) > max_windows)
  regions <- regions[seq_len(max_windows), , drop = FALSE]
log_msg(sprintf("Starting native STAAR window chunk %d [mode=%s, run %s]: chr%s [%d, %d], %d windows",
                chunk_id, spa_mode, run_name, chr, chunk$chunk_start, chunk$chunk_end, nrow(regions)))
log_msg(sprintf("Config snapshot from %s; maf=%s rv_cutoff=%d use_SPA(null)=%s",
                config$prepared_at, config$filter_spec$rare_maf_cutoff,
                null_or(config$staar_rv_num_cutoff, 2L), null_or(obj_nullmodel$use_SPA, "NULL")))

# ---- Open the chr GDS ----
gds_path <- file.path(config$gds_dir, gsub("\\{chr\\}", chr, config$gds_pattern))
stopifnot(file.exists(gds_path))
genofile <- seqOpen(gds_path)
on.exit(seqClose(genofile), add = TRUE)

# ---- Positional extraction of the native row (layout asserted; see header) ----
.curated_row <- function(res_row, label, mode) {
  nms <- colnames(res_row); nc <- length(res_row)
  stopifnot(identical(nms[1:5], c("Chr", "Start Loc", "End Loc", "#SNV", "cMAC")))
  v <- suppressWarnings(as.numeric(res_row))
  base <- data.frame(label = label, chr = as.character(v[1L]),
                     start = as.integer(v[2L]), end = as.integer(v[3L]),
                     stringsAsFactors = FALSE)
  if (identical(mode, "nospa")) {
    stopifnot(identical(nms[(nc - 1L):nc], c("ACAT-O", "STAAR-O")),
              (nc - 7L) %% 6L == 0L)
    B <- (nc - 7L) %/% 6L                    # family block length (last = family CCT)
    fam <- function(m) v[5L + m * B]
    base$n_SNV_native       <- as.integer(v[4L])
    base$cMAC_native        <- v[5L]
    base$SKAT_1_25_native   <- fam(1L)
    base$SKAT_1_1_native    <- fam(2L)
    base$Burden_1_25_native <- fam(3L)
    base$Burden_1_1_native  <- fam(4L)
    base$ACATV_1_25_native  <- fam(5L)
    base$ACATV_1_1_native   <- fam(6L)
    base$ACAT_O_native      <- v[nc - 1L]
    base$STAAR_O_native     <- v[nc]
  } else {
    stopifnot(identical(nms[nc], "STAAR-B"), (nc - 6L) %% 2L == 0L)
    B <- (nc - 6L) %/% 2L
    base$n_SNV_spa       <- as.integer(v[4L])
    base$cMAC_spa        <- v[5L]
    base$Burden_1_25_spa <- v[5L + B]
    base$Burden_1_1_spa  <- v[5L + 2L * B]
    base$STAAR_B_native  <- v[nc]
  }
  base
}
.schema_cols <- function(mode) {
  if (identical(mode, "nospa"))
    c("label","chr","start","end","n_SNV_native","cMAC_native",
      "SKAT_1_25_native","SKAT_1_1_native","Burden_1_25_native","Burden_1_1_native",
      "ACATV_1_25_native","ACATV_1_1_native","ACAT_O_native","STAAR_O_native")
  else c("label","chr","start","end","n_SNV_spa","cMAC_spa",
         "Burden_1_25_spa","Burden_1_1_spa","STAAR_B_native")
}

# ---- Loop windows: fully native Sliding_Window(type = "single") per window ----
rows <- vector("list", nrow(regions))
n_ok <- n_skip <- n_untested <- n_error <- 0L
for (k in seq_len(nrow(regions))) {
  w <- regions[k, ]
  out <- tryCatch(
    STAARpipeline::Sliding_Window(
      chr = as.numeric(chr), start_loc = w$start, end_loc = w$end, type = "single",
      genofile = genofile, obj_nullmodel = obj_nullmodel,
      rare_maf_cutoff = config$filter_spec$rare_maf_cutoff,
      rv_num_cutoff   = null_or(config$staar_rv_num_cutoff, 2L),
      QC_label        = null_or(config$filter_spec$qc_label, "annotation/filter"),
      variant_type    = config$filter_spec$variant_type,
      geno_missing_imputation = null_or(config$geno_missing_imputation, "mean"),
      Annotation_dir  = null_or(config$Annotation_dir,
                                "annotation/info/FunctionalAnnotation"),
      Annotation_name_catalog = staar_catalog,
      Use_annotation_weights  = use_anno, Annotation_name = Annotation_name,
      SPA_p_filter = FALSE, silent = TRUE),
    error = function(e) e)
  if (inherits(out, "error")) {
    msg <- conditionMessage(out)
    if (grepl("less than 2", msg, fixed = TRUE)) {
      # Routine skip: fewer than 2 PASS variants in the window.
      n_skip <- n_skip + 1L
    } else if (grepl("object 'results' not found", msg, fixed = TRUE)) {
      # Sliding_Window_Single's failure path: STAAR() itself errored inside
      # (fewer than 2 testable variants after its internal MAF re-filter, or a
      # numerical failure), so its `results` was never defined. Counted as
      # untested-by-native; the first few are logged for diagnosis.
      n_untested <- n_untested + 1L
      if (n_untested <= 5L)
        log_msg(sprintf("  window %s untested by native STAAR (internal failure)", w$label))
    } else {
      n_error <- n_error + 1L
      log_msg(sprintf("  window %s ERROR: %s", w$label, msg))
    }
    next
  }
  if (is.null(out) || nrow(out) == 0L) { n_skip <- n_skip + 1L; next }
  rows[[k]] <- .curated_row(out[1L, , drop = FALSE], w$label, spa_mode)
  n_ok <- n_ok + 1L
  if (k %% 100L == 0L)
    log_msg(sprintf("  ... %d/%d windows (ok=%d skip=%d untested=%d err=%d, %.1fs)",
                    k, nrow(regions), n_ok, n_skip, n_untested, n_error,
                    as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
}
rows <- rows[!vapply(rows, is.null, logical(1L))]
curated <- if (length(rows)) do.call(rbind, rows) else
  as.data.frame(setNames(rep(list(character(0L)), length(.schema_cols(spa_mode))),
                         .schema_cols(spa_mode)), stringsAsFactors = FALSE)

# ---- Persist ----
out_path <- file.path(results_staar_dir, sprintf("staar_window_%04d.csv", chunk_id))
write.csv(curated, out_path, row.names = FALSE)
log_msg(sprintf("Wrote %s (%d window rows)", out_path, nrow(curated)))
log_msg(sprintf("Done native STAAR window chunk %d [mode=%s, annotation_weights=%s]. windows=%d, n_ok=%d, n_skip=%d, n_untested=%d, n_error=%d, wall=%.1fs",
                chunk_id, spa_mode, use_anno, nrow(regions), n_ok, n_skip,
                n_untested, n_error,
                as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
