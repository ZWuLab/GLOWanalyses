# ============================================================================
# _dataprep_lib.R - shared helpers for the 00-data-prep templates
#
# Sourced (not run) by plink-to-gds.R / annotate-favor.R / compute-pcs.R /
# assemble-pheno-covar.R. Encodes the cohort-data LINEAGE convention: data-prep
# outputs live in a user-designated `data_root` as a sibling-suffix tree, NOT
# under runs/.
#
#   <data_root>/<base_name>_plink/      input PLINK (chr{chr}...)
#   <data_root>/<base_name>_gds/        plink_to_gds      -> chr{chr}...gds
#   <data_root>/<base_name>_gds_favor/  annotate_favor    -> {match}/{gds,csv}/
#   <data_root>/<base_name>_pcs/        compute_pcs_gds   -> pcs.rds/.csv
#   <data_root>/<base_name>_pheno/      assemble_pheno_covar -> <run>_pheno_covar.rds
#
# Each step writes a README provenance header + a resolved-config snapshot +
# a logs/ dir into its output directory, so every derived dataset is
# self-documenting and reproducible (codifying the existing processed/als/
# convention). File NAMES (e.g. chr{chr}_hg38_favor.gds) stay in the config as
# patterns - only the directory layout is standardized here.
#
# SHARED STORAGE (not per-run born-shared): data-prep outputs are cohort DATA
# ASSETS in the user's `data_root` lineage tree (often OUTSIDE the repo), not
# run-org outputs/ leaves, so the `symlinked_shared_root` per-run born-sharing
# convention does NOT apply here. Share these at the `data_root` level instead
# (place that root on a shared mount).

null_or <- function(a, b) if (is.null(a)) b else a

# Parse `--config <path>` (+ pass through any extra `--flag value` pairs as a
# named character list). Returns list(config = <path>, opts = <named list>).
parse_dataprep_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  config_path <- NA_character_; opts <- list(); i <- 1L
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--config" && i < length(args)) { config_path <- args[i + 1L]; i <- i + 2L }
    else if (startsWith(a, "--") && i < length(args)) {
      opts[[sub("^--", "", a)]] <- args[i + 1L]; i <- i + 2L
    } else stop("Unknown / dangling argument: ", a)
  }
  if (is.na(config_path)) stop("Missing required --config <run>/config.R")
  # Relative paths resolve from the GLOWanalyses directory (the run-from convention);
  # name the path and the CWD so a typo'd config is obvious in a batch log.
  if (!file.exists(config_path))
    stop(sprintf("Run config not found: %s  (working directory: %s)", config_path, getwd()))
  list(config = config_path, opts = opts)
}

# NOTE on config sourcing: templates source the run config at their TOP LEVEL with
#   source(pa$config, local = TRUE); cfg <- environment()
# (the proven 03-snv-set idiom) so the config and its own source(base, local=TRUE)
# land in the same frame and base functions resolve. Read fields with
#   g0 <- function(nm, d = NULL) get0(nm, envir = cfg, ifnotfound = d, inherits = FALSE)

# Standard sibling-suffix directories for a cohort lineage tree.
resolve_cohort_paths <- function(data_root, base_name, match_method = NULL) {
  stopifnot(nzchar(data_root), nzchar(base_name))
  sib <- function(suffix) file.path(data_root, paste0(base_name, suffix))
  paths <- list(
    data_root  = data_root,
    base_name  = base_name,
    plink_dir  = sib("_plink"),
    gds_dir    = sib("_gds"),
    favor_dir  = sib("_gds_favor"),
    pcs_dir    = sib("_pcs"),
    pheno_dir  = sib("_pheno"))
  if (!is.null(match_method)) {
    paths$favor_gds_dir <- file.path(paths$favor_dir, match_method, "gds")
    paths$favor_csv_dir <- file.path(paths$favor_dir, match_method, "csv")
  }
  paths
}

# Expand a per-chromosome file pattern (with a literal "{chr}" token).
chr_path <- function(dir, pattern, chr) file.path(dir, gsub("{chr}", chr, pattern, fixed = TRUE))

# Write a self-documenting provenance record into a data-prep output dir:
#   <dir>/README.md            (header matching the processed/als convention)
#   <dir>/provenance/config_snapshot.rds + sessionInfo.txt
#   <dir>/logs/                (created)
# `unit`: when a step runs as one task of a per-chromosome array (`--chr N`),
# pass e.g. "chr22" so the task writes provenance/config_snapshot_chr22.rds +
# sessionInfo_chr22.txt (one record per task) instead of every task overwriting
# the single serial-run snapshot. The README is shared: it is identical across
# tasks and written atomically so concurrent tasks cannot interleave it.
write_dataprep_provenance <- function(dir, step, script, source_alias = NA, source_path = NA,
                                      config_snapshot = NULL, notes = NULL, unit = NULL) {
  dir.create(file.path(dir, "provenance"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(dir, "logs"), recursive = TRUE, showWarnings = FALSE)
  suffix <- if (is.null(unit)) "" else paste0("_", unit)
  # ignore.stderr: on compute nodes without git on PATH the subshell's
  # "sh: 1: git: not found" would otherwise leak into the SLURM .err log.
  git_sha <- tryCatch(suppressWarnings(system("git rev-parse --short HEAD",
                      intern = TRUE, ignore.stderr = TRUE)),
                      error = function(e) NA_character_)
  if (length(git_sha) != 1L || is.na(git_sha) || !nzchar(git_sha)) git_sha <- "unknown"
  readme <- c(
    sprintf("# %s", basename(dir)), "",
    sprintf("**Step**: %s", step),
    sprintf("**Source alias**: %s", source_alias),
    sprintf("**Source path**: %s", source_path),
    sprintf("**Script**: %s", script),
    sprintf("**Git SHA**: %s", git_sha), "",
    if (!is.null(notes)) c("## Notes", "", notes, "") else NULL,
    "Generated by a GLOWanalyses 00-data-prep template; see provenance/ for the resolved config",
    "(one config_snapshot_<unit>.rds + sessionInfo_<unit>.txt per task when run as a per-chromosome array).")
  # Atomic README write: temp file in the same dir, then rename (a POSIX rename is
  # atomic), so array tasks finishing together cannot interleave their writes.
  tmp <- tempfile(pattern = "README_", tmpdir = dir, fileext = ".md.tmp")
  writeLines(readme, tmp)
  file.rename(tmp, file.path(dir, "README.md"))
  if (!is.null(config_snapshot))
    saveRDS(config_snapshot, file.path(dir, "provenance", paste0("config_snapshot", suffix, ".rds")))
  writeLines(capture.output(sessionInfo()),
             file.path(dir, "provenance", paste0("sessionInfo", suffix, ".txt")))
  invisible(file.path(dir, "README.md"))
}
