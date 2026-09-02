#!/usr/bin/env Rscript
# ============================================================================
# 03-snv-set/run-staar-coding.R - Native STAARpipeline coding engine
#
# Faithful adaptation of legacy-materials/STAARpipeline-Tutorial/
# STAARpipeline_Gene_Centric_Coding.r for the ALS essentialdb GDS. Gene-chunk
# SLURM array: `arrayid` -> (chr, gene-group) via the tutorial's
# ceiling(table(genes_info[,2]) / gene_num_in_array) mapping. Run ONCE PER MODE
# (`--spa {on,off}`); the mode selects which Stage-1 STAAR null model to load
# (Gene_Centric_Coding auto-detects use_SPA) and which output dir to write:
#   --spa on  -> 01-staar_null_model_spa.rds   -> results-staar-spa/   -> STAAR_B_native
#   --spa off -> 01-staar_null_model_nospa.rds -> results-staar-nospa/ -> STAAR_O_native
#
# Each Gene_Centric_Coding(category = "all_categories_incl_ptv", ...) call reads a
# gene once and returns ALL 7 categories (no per-category re-extraction). The
# essentialdb Annotation_name_catalog (Stage-1) maps the masking + 11 weighting
# nodes; native annotation imputation is internal (no external medians needed).
#
# Outputs (under <output_dir>/):
#   results-staar-<mode>/staar_coding_<arrayid>.csv   curated flat per-(gene, category) table
#       (mode-keyed: nospa -> STAAR_O_native + ACAT_O + the SKAT/Burden/ACAT-V grid;
#       spa -> STAAR_B_native + burden), genome-rbind'd by Stage 3. A 0-row schema frame
#       when the task's genes returned nothing.
#   staar_detail/staar_coding_full_<mode>_<arrayid>.csv  OPT-IN (config$staar_write_subtests):
#       the full annotation-weighted STAAR grid (every SKAT/Burden/ACAT-V x annotation sub-test).
#   logs/log_staar_<mode>_<arrayid>.txt               per-task summary + errors
#
# Usage (from the GLOWanalyses directory):
#   conda activate r_env
#   Rscript 03-snv-set/run-staar-coding.R \
#       --arrayid <A> \
#       --config runs/example/config.R \
#       --spa on

suppressMessages({
  library(gdsfmt); library(SeqArray); library(SeqVarTools)
  library(STAAR); library(STAARpipeline)
})

null_or <- function(a, b) if (is.null(a)) b else a
t_start <- Sys.time()

# ---- Args ----
args <- commandArgs(trailingOnly = TRUE)
arrayid <- NULL; config_path <- NA_character_; spa_on <- NA
i <- 1L
while (i <= length(args)) {
  if (args[i] == "--arrayid" && i < length(args)) { arrayid <- as.integer(args[i + 1L]); i <- i + 2L }
  else if (args[i] == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
  else if (args[i] == "--spa" && i < length(args)) {
    v <- tolower(args[i + 1L])
    spa_on <- if (v %in% c("on", "true", "1")) {
      TRUE
    } else if (v %in% c("off", "false", "0")) {
      FALSE
    } else {
      stop("--spa must be on|off (got '", args[i + 1L], "')")
    }
    i <- i + 2L
  } else stop("Unknown argument: ", args[i])
}
stopifnot(!is.null(arrayid), arrayid >= 1L)
if (is.na(config_path)) stop("Missing required --config <run>/config.R")
if (is.na(spa_on)) stop("Missing required --spa {on,off}")
stopifnot(file.exists(config_path))
spa_mode <- if (spa_on) "spa" else "nospa"

# ---- Run identity: output_dir DERIVED from --config ----
output_dir <- file.path(dirname(config_path), "outputs")
run_name   <- basename(dirname(config_path))
shared_dir         <- file.path(output_dir, "shared")
results_staar_dir  <- file.path(output_dir, paste0("results-staar-", spa_mode))
logs_dir           <- file.path(output_dir, "logs")
for (d in c(results_staar_dir, logs_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Load Stage-1 artifacts ----
config <- readRDS(file.path(shared_dir, "01-config_used.rds"))
if (!isTRUE(config$staar_enabled)) stop("Stage 1 says staar_enabled = FALSE; nothing to run.")
null_path <- file.path(shared_dir, if (spa_on) "01-staar_null_model_spa.rds"
                                   else        "01-staar_null_model_nospa.rds")
catalog_path <- file.path(shared_dir, "01-staar_annotation_catalog.rds")
for (p in c(null_path, catalog_path)) {
  if (!file.exists(p)) stop("Missing Stage-1 artifact: ", p, " (run 01-prepare-shared.R).")
}
obj_nullmodel <- readRDS(null_path)        # use_SPA auto-detected by Gene_Centric_Coding
staar_catalog <- readRDS(catalog_path)
Annotation_name <- null_or(config$staar_anno_features, config$pi_features)
gene_num_in_array <- config$gene_num_in_array

# ---- Logging ----
log_path <- file.path(logs_dir, sprintf("log_staar_%s_%04d.txt", spa_mode, arrayid))
log_con  <- file(log_path, open = "wt")
on.exit(try(close(log_con), silent = TRUE), add = TRUE, after = FALSE)
log_msg <- function(msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n", sep = ""); writeLines(line, log_con); flush(log_con)
}

# ---- arrayid -> (chr, gene-group)  [tutorial logic; genes_info is the STAARpipeline global] ----
group.num.allchr <- ceiling(table(genes_info[, 2]) / gene_num_in_array)
n_array_total    <- sum(group.num.allchr)
if (arrayid > n_array_total) {
  stop(sprintf("arrayid %d exceeds the gene-chunk count (%d at %d genes/array).",
               arrayid, n_array_total, gene_num_in_array))
}
chr     <- which.max(arrayid <= cumsum(group.num.allchr))
groupid <- if (chr == 1L) arrayid else arrayid - cumsum(group.num.allchr)[chr - 1L]
genes_chr <- genes_info[genes_info[, 2] == chr, , drop = FALSE]
sub_seq_num <- nrow(genes_chr)
ids <- if (groupid < group.num.allchr[chr]) {
  ((groupid - 1L) * gene_num_in_array + 1L):(groupid * gene_num_in_array)
} else {
  ((groupid - 1L) * gene_num_in_array + 1L):sub_seq_num
}

log_msg(sprintf("Starting STAAR coding [mode=%s, run %s]: arrayid=%d -> chr%d group %d (%d genes)",
                spa_mode, run_name, arrayid, chr, groupid, length(ids)))
log_msg(sprintf("Config snapshot from %s; maf=%s rv_cutoff=%d use_SPA(null)=%s",
                config$prepared_at, config$filter_spec$rare_maf_cutoff,
                config$staar_rv_num_cutoff, null_or(obj_nullmodel$use_SPA, "NULL")))

# ---- Open the chr GDS ----
gds_path <- file.path(config$gds_dir, gsub("\\{chr\\}", chr, config$gds_pattern))
stopifnot(file.exists(gds_path))
genofile <- seqOpen(gds_path)
on.exit(seqClose(genofile), add = TRUE)

# ---- Loop genes; Gene_Centric_Coding all 7 categories per gene ----
# SPA_p_filter = FALSE is REQUIRED (Phase-0): the default TRUE forwards a
# marginal-p pre-filter into STAAR_Binary_SPA that degenerates the genotype matrix
# on small variant sets -> "Not a matrix.". FALSE runs the complete set test (the
# pre-filter is only a speed optimization, negligible on this sparse cohort).
res <- vector("list", length(ids)); names_acc <- character(0L)
n_ok <- n_empty <- n_error <- 0L
for (kk in ids) {
  gene_name <- as.character(genes_chr[kk, 1])
  out <- tryCatch(
    Gene_Centric_Coding(
      chr = chr, gene_name = gene_name, category = "all_categories_incl_ptv",
      genofile = genofile, obj_nullmodel = obj_nullmodel,
      rare_maf_cutoff = config$filter_spec$rare_maf_cutoff,
      rv_num_cutoff   = config$staar_rv_num_cutoff,
      QC_label        = null_or(config$filter_spec$qc_label, "annotation/filter"),
      variant_type    = config$filter_spec$variant_type,
      geno_missing_imputation = config$geno_missing_imputation,
      Annotation_dir  = config$Annotation_dir,
      Annotation_name_catalog = staar_catalog,
      Use_annotation_weights  = TRUE, Annotation_name = Annotation_name,
      SPA_p_filter = FALSE),
    error = function(e) { log_msg(sprintf("  gene %s ERROR: %s", gene_name, conditionMessage(e))); structure("error", class = "staar_err") })
  if (inherits(out, "staar_err")) { n_error <- n_error + 1L; next }
  # Keep only genes that returned >=1 non-empty category.
  non_empty <- !is.null(out) && length(out) > 0L &&
    any(vapply(out, function(x) !is.null(x) && nrow(as.data.frame(x)) > 0L, logical(1)))
  if (non_empty) {
    res[[gene_name]] <- out; names_acc <- c(names_acc, gene_name); n_ok <- n_ok + 1L
  } else n_empty <- n_empty + 1L
}
res <- res[names_acc]    # drop unfilled slots; keep only genes with results

# ---- Flatten (the per-task output is now a flat CSV, human-reviewable) ----
# Pull a raw STAAR column by exact name with NA fallback (the per-mode return omits
# some columns, and the missense-family categories carry extra "-Disruptive" cols).
.pull_num <- function(df, col) if (col %in% names(df)) suppressWarnings(as.numeric(df[[col]][1L])) else NA_real_

# Curated flat table (always written): the headline omnibi + base grid that Stage 3
# aggregates, mode-keyed (_native / _spa). Returns a 0-row frame with the mode schema
# when the task found nothing, so Stage 3 can rbind uniformly.
empty_df <- function(cols) as.data.frame(setNames(rep(list(character(0L)), length(cols)), cols),
                                         stringsAsFactors = FALSE)
flatten_curated <- function(res, mode) {
  rows <- list()
  for (gene in names(res)) {
    glist <- res[[gene]]; if (is.null(glist)) next
    for (category in names(glist)) {
      df <- glist[[category]]; if (is.null(df)) next
      df <- as.data.frame(df, check.names = FALSE); if (nrow(df) == 0L) next
      base <- data.frame(gene = gene, chr = as.character(df[["Chr"]][1L]), category = category,
                         stringsAsFactors = FALSE)
      if (identical(mode, "nospa")) {
        base$n_SNV_native       <- as.integer(df[["#SNV"]][1L])
        base$cMAC_native        <- .pull_num(df, "cMAC")
        base$STAAR_O_native     <- .pull_num(df, "STAAR-O")
        base$ACAT_O_native      <- .pull_num(df, "ACAT-O")
        base$SKAT_1_25_native   <- .pull_num(df, "SKAT(1,25)")
        base$SKAT_1_1_native    <- .pull_num(df, "SKAT(1,1)")
        base$Burden_1_25_native <- .pull_num(df, "Burden(1,25)")
        base$Burden_1_1_native  <- .pull_num(df, "Burden(1,1)")
        base$ACATV_1_25_native  <- .pull_num(df, "ACAT-V(1,25)")
        base$ACATV_1_1_native   <- .pull_num(df, "ACAT-V(1,1)")
      } else {  # spa: burden-only -> STAAR-B
        base$n_SNV_spa       <- as.integer(df[["#SNV"]][1L])
        base$cMAC_spa        <- .pull_num(df, "cMAC")
        base$STAAR_B_native  <- .pull_num(df, "STAAR-B")
        base$Burden_1_25_spa <- .pull_num(df, "Burden(1,25)")
        base$Burden_1_1_spa  <- .pull_num(df, "Burden(1,1)")
      }
      rows[[length(rows) + 1L]] <- base
    }
  }
  if (length(rows) == 0L) {
    cols <- if (identical(mode, "nospa"))
      c("gene","chr","category","n_SNV_native","cMAC_native","STAAR_O_native","ACAT_O_native",
        "SKAT_1_25_native","SKAT_1_1_native","Burden_1_25_native","Burden_1_1_native",
        "ACATV_1_25_native","ACATV_1_1_native")
    else c("gene","chr","category","n_SNV_spa","cMAC_spa","STAAR_B_native","Burden_1_25_spa","Burden_1_1_spa")
    return(empty_df(cols))
  }
  do.call(rbind, rows)
}

# Full annotation-weighted grid (opt-in): every raw STAAR column (SKAT/Burden/ACAT-V x
# annotation sub-test), union-aligned across categories (the "-Disruptive" cols are
# absent from synonymous etc.) with NA fill; widest schema leads the column order.
# STAAR returns each cell as a length-1 LIST element (every column is a list), which
# write.csv cannot encode, so each column is first coerced to its atomic scalar type
# (length-preserving; NULL/empty -> NA).
.atomize_col <- function(col) {
  out <- vector("list", length(col))
  for (i in seq_along(col)) out[[i]] <- if (length(col[[i]])) col[[i]][[1L]] else NA
  unlist(out, use.names = FALSE)
}
.atomize_df <- function(df) as.data.frame(lapply(df, .atomize_col),
                                          check.names = FALSE, stringsAsFactors = FALSE)
flatten_full <- function(res) {
  dfs <- list()
  for (gene in names(res)) {
    glist <- res[[gene]]; if (is.null(glist)) next
    for (category in names(glist)) {
      df <- glist[[category]]; if (is.null(df)) next
      df <- as.data.frame(df, check.names = FALSE); if (nrow(df) == 0L) next
      dfs[[length(dfs) + 1L]] <- .atomize_df(df)
    }
  }
  if (length(dfs) == 0L) return(NULL)
  dfs <- dfs[order(vapply(dfs, ncol, integer(1)), decreasing = TRUE)]
  all_cols <- Reduce(union, lapply(dfs, names))
  aligned <- lapply(dfs, function(d) { miss <- setdiff(all_cols, names(d))
    if (length(miss)) d[miss] <- NA; d[all_cols] })
  do.call(rbind, aligned)
}

# ---- Persist (flat CSV; curated always, full annotation grid opt-in) ----
curated  <- flatten_curated(res, spa_mode)
out_path <- file.path(results_staar_dir, sprintf("staar_coding_%04d.csv", arrayid))
write.csv(curated, out_path, row.names = FALSE)
log_msg(sprintf("Wrote %s (%d gene-category rows; %d genes with results)",
                out_path, nrow(curated), length(res)))

if (isTRUE(config$staar_write_subtests)) {
  full <- flatten_full(res)
  if (!is.null(full)) {
    detail_dir <- file.path(output_dir, "staar_detail")
    dir.create(detail_dir, recursive = TRUE, showWarnings = FALSE)
    full_path <- file.path(detail_dir, sprintf("staar_coding_full_%s_%04d.csv", spa_mode, arrayid))
    write.csv(full, full_path, row.names = FALSE)
    log_msg(sprintf("Wrote %s (%d rows x %d cols; full annotation-weighted grid)",
                    full_path, nrow(full), ncol(full)))
  }
}
log_msg(sprintf("Done STAAR coding arrayid=%d [mode=%s]. genes=%d, n_ok=%d, n_empty=%d, n_error=%d, wall=%.1fs",
                arrayid, spa_mode, length(ids), n_ok, n_empty, n_error,
                as.numeric(difftime(Sys.time(), t_start, units = "secs"))))
