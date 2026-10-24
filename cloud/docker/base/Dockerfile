# syntax=docker/dockerfile:1
# ==============================================================================
# File: cloud/docker/base/Dockerfile
# Purpose: the GLOW base image, glow-base. R 4.3.3 from rocker/r-ver with the dated CRAN
#   snapshot that image configures (p3m.dev, 2024-04-23), Bioconductor 3.18, and the STAAR
#   family at the commits installed on the authors' HPC, so the image stays on the environment
#   in which the GLOW results were computed. It is rebuilt only when this list changes and is
#   tagged by date (GLOW_BASE_TAG in ../versions.env). The GLOW packages and the GLOWanalyses
#   scripts are added by ../Dockerfile. Pins: STAAR 0.9.8, MultiSTAAR 0.9.7.1, SCANG 1.0.4,
#   STAARpipeline 0.9.8 (the last 0.9.8 commit; identical DESCRIPTION and src/ to the copy the
#   authors validated against). Rcpp and RcppArmadillo come from the same snapshot as every
#   other CRAN package (RcppArmadillo 0.12.8.2.1): GFisher 0.3.2 compiles against it, 0.3.1 did not.
#   Authors: Zheyang Wu (Worcester Polytechnic Institute), for ZWuLab.
#   Derived from the Terra workflows written by Alberto Brusati (Landers lab, UMass Chan Medical
#   School), 2026: his Dockerfile for abrusati/glow-r:2.0 set the shape of this recipe (archives
#   at pinned commits, no GitHub API calls, a smoke test that loads everything).
FROM rocker/r-ver:4.3.3

LABEL org.opencontainers.image.source="https://github.com/ZWuLab/GLOWanalyses" \
      org.opencontainers.image.description="GLOW base image: R 4.3.3, Bioconductor 3.18, the STAAR family" \
      org.opencontainers.image.licenses="GPL-3.0"

# The pinned commits (GitHub archives; a plain download, no API call).
ARG STAAR_SHA=9db9dd504905b9f469146f670e5f6dbe3e08d01a
ARG MULTISTAAR_SHA=c372e135d88d5537c43af2d0f3e935f47cafd11c
ARG SCANG_SHA=89f16cf21b04132cb9cdce2940c03a119411941a
ARG STAARPIPELINE_SHA=7e868f23c85ba0713dd00c8512ad23cca7fbc8f9
ARG BIOC_VERSION=3.18

# System libraries: curl, ssl and xml2 for R's downloads and the XML parsers that rtracklayer
# and GenomicFeatures pull in; zlib, bz2 and lzma for gdsfmt, SeqArray and Rhtslib; git so that
# GLOWpipeline's repository detection has the binary it calls (it falls back cleanly without it).
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev libbz2-dev liblzma-dev \
    libpng-dev git ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# CRAN packages, as binaries from the snapshot rocker configured (printed first for the log).
# arrow reads the FAVOR 2.0 Parquet files; its Parquet support is asserted.
RUN Rscript -e ' \
  cat("CRAN:", getOption("repos")[["CRAN"]], "\n"); \
  options(Ncpus = max(1L, parallel::detectCores())); \
  install.packages(c("BiocManager", "Rcpp", "RcppArmadillo", "mvtnorm", "glmnet", "Matrix", \
                     "SPAtest", "digest", "data.table", "dplyr", "fst", "readxl", "MASS", \
                     "kinship2", "arrow")); \
  for (p in c("Rcpp", "RcppArmadillo", "mvtnorm", "glmnet", "Matrix", "SPAtest", "digest", \
              "data.table", "dplyr", "fst", "readxl", "kinship2", "arrow")) \
    if (!requireNamespace(p, quietly = TRUE)) stop("CRAN package missing: ", p); \
  stopifnot(arrow::arrow_info()$capabilities[["parquet"]]); \
  cat("arrow", as.character(packageVersion("arrow")), "with Parquet\n")'

# Bioconductor 3.18 (the release for R 4.3). GMMAT is a CRAN package whose dependencies are
# Bioconductor packages, so it is installed here, through BiocManager, which reads both.
RUN MAKEFLAGS=-j4 Rscript -e ' \
  options(Ncpus = max(1L, parallel::detectCores())); \
  BiocManager::install(version = Sys.getenv("BIOC_VERSION"), ask = FALSE, update = FALSE); \
  pkgs <- c("gdsfmt", "SeqArray", "SeqVarTools", "SNPRelate", "GenomicRanges", "GenomicFeatures", \
            "TxDb.Hsapiens.UCSC.hg38.knownGene", "GENESIS", "GMMAT"); \
  BiocManager::install(pkgs, version = Sys.getenv("BIOC_VERSION"), ask = FALSE, update = FALSE); \
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) stop("Bioconductor package missing: ", p); \
  cat("Bioconductor", as.character(BiocManager::version()), "\n")'

# The STAAR family from GitHub archives at the pinned commits, in dependency order
# (MultiSTAAR and SCANG link against STAAR's headers). Installed from source, unmodified.
RUN MAKEFLAGS=-j4 Rscript -e ' \
  inst <- function(url, pkg) { \
    tgz <- tempfile(fileext = ".tar.gz"); download.file(url, tgz, quiet = TRUE); \
    ex <- tempfile(); untar(tgz, exdir = ex); \
    install.packages(list.dirs(ex, recursive = FALSE)[1], repos = NULL, type = "source"); \
    if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " did not install"); \
    cat(pkg, as.character(packageVersion(pkg)), "from", url, "\n") }; \
  inst(sprintf("https://github.com/xihaoli/STAAR/archive/%s.tar.gz", Sys.getenv("STAAR_SHA")), "STAAR"); \
  inst(sprintf("https://github.com/xihaoli/MultiSTAAR/archive/%s.tar.gz", Sys.getenv("MULTISTAAR_SHA")), "MultiSTAAR"); \
  inst(sprintf("https://github.com/zilinli1988/SCANG/archive/%s.tar.gz", Sys.getenv("SCANG_SHA")), "SCANG"); \
  inst(sprintf("https://github.com/xihaoli/STAARpipeline/archive/%s.tar.gz", Sys.getenv("STAARPIPELINE_SHA")), "STAARpipeline")'

# Smoke test: everything loads; the package manifest is kept in the image for provenance.
RUN Rscript -e ' \
  for (p in c("SeqArray", "gdsfmt", "SeqVarTools", "SNPRelate", "GenomicFeatures", \
              "TxDb.Hsapiens.UCSC.hg38.knownGene", "GENESIS", "GMMAT", "STAAR", "MultiSTAAR", \
              "SCANG", "STAARpipeline", "arrow", "glmnet", "SPAtest", "data.table")) \
    suppressPackageStartupMessages(library(p, character.only = TRUE)); \
  ip <- installed.packages()[, c("Package", "Version")]; \
  write.table(ip, "/opt/glow-base-manifest.tsv", sep = "\t", quote = FALSE, row.names = FALSE); \
  cat(sprintf("glow-base: %s | %s | RcppArmadillo %s | %d packages\n", R.version.string, \
      paste("Bioconductor", as.character(BiocManager::version())), \
      as.character(packageVersion("RcppArmadillo")), nrow(ip)))'
