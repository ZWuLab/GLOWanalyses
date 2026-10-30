# syntax=docker/dockerfile:1
# ==============================================================================
# File: cloud/docker/Dockerfile
# Purpose: the GLOW image, glowanalyses, built on glow-base (base/Dockerfile): GFisher, GLOWr
#   and GLOWpipeline installed from their release-tag archives at the versions named in
#   versions.env, and this GLOWanalyses tree copied to /opt/GLOWanalyses, the working directory,
#   so that a WDL task runs the templates exactly as the README says. Build it with the
#   GLOWanalyses root as the context and this file through -f; the recipe in
#   .github/workflows/docker-image.yaml does so at every release tag and tags the image with it.
#   Authors: Zheyang Wu (Worcester Polytechnic Institute), for ZWuLab.
#   Derived from the Terra workflows written by Alberto Brusati (Landers lab, UMass Chan Medical
#   School), 2026: the archive-install pattern with version asserts is his.
ARG BASE_IMAGE=ghcr.io/zwulab/glow-base:2026-10-22.4
FROM ${BASE_IMAGE}

ARG GFISHER_TAG=v0.3.2
ARG GLOWR_VERSION=0.3.4
ARG GLOWPIPELINE_VERSION=0.2.2
ARG GLOWANALYSES_REF=unknown

LABEL org.opencontainers.image.source="https://github.com/ZWuLab/GLOWanalyses" \
      org.opencontainers.image.description="GLOW: GLOWanalyses templates with GFisher, GLOWr and GLOWpipeline" \
      org.opencontainers.image.licenses="GPL-3.0" \
      org.opencontainers.image.version="${GLOWANALYSES_REF}"

# The three GLOW packages from the tagged archives (a plain download, no API call), in
# dependency order, each asserted at the pinned version. GFisher installs unmodified.
RUN MAKEFLAGS=-j4 Rscript -e ' \
  inst <- function(url, pkg, version) { \
    tgz <- tempfile(fileext = ".tar.gz"); download.file(url, tgz, quiet = TRUE); \
    ex <- tempfile(); untar(tgz, exdir = ex); \
    install.packages(list.dirs(ex, recursive = FALSE)[1], repos = NULL, type = "source"); \
    if (!requireNamespace(pkg, quietly = TRUE)) stop(pkg, " did not install"); \
    v <- as.character(packageVersion(pkg)); \
    if (v != version) stop(sprintf("%s is %s, expected %s", pkg, v, version)); \
    cat(pkg, v, "installed\n") }; \
  gf <- Sys.getenv("GFISHER_TAG"); \
  inst(sprintf("https://github.com/ZWuLab/GFisher/archive/refs/tags/%s.tar.gz", gf), "GFisher", sub("^v", "", gf)); \
  inst(sprintf("https://github.com/ZWuLab/GLOWr/archive/refs/tags/v%s.tar.gz", Sys.getenv("GLOWR_VERSION")), "GLOWr", Sys.getenv("GLOWR_VERSION")); \
  inst(sprintf("https://github.com/ZWuLab/GLOWpipeline/archive/refs/tags/v%s.tar.gz", Sys.getenv("GLOWPIPELINE_VERSION")), "GLOWpipeline", Sys.getenv("GLOWPIPELINE_VERSION"))'

# The GLOWanalyses tree (the build context, minus .dockerignore) at its released tag.
COPY . /opt/GLOWanalyses
WORKDIR /opt/GLOWanalyses

# One BLAS thread per process: the templates parallelize with forked workers (n_cores), so a
# threaded BLAS would oversubscribe the machine. Rscript is the one rocker ships.
ENV OPENBLAS_NUM_THREADS=1 \
    GLOW_RSCRIPT=/usr/local/bin/Rscript

# Smoke test: the GLOW stack loads, the tree is in place, and the manifest records what is here.
RUN Rscript -e ' \
  for (p in c("GFisher", "GLOWr", "GLOWpipeline", "SeqArray", "STAAR")) \
    suppressPackageStartupMessages(library(p, character.only = TRUE)); \
  stopifnot(file.exists("run-example.sh"), file.exists("03-snv-set/prepare.R")); \
  src <- if (file.exists(".release-source")) paste(readLines(".release-source"), collapse = " | ") else "(no .release-source)"; \
  lines <- c(sprintf("GLOWanalyses ref: %s", Sys.getenv("GLOWANALYSES_REF")), sprintf("release-source: %s", src), \
             sprintf("%s %s", c("GFisher", "GLOWr", "GLOWpipeline", "STAAR", "STAARpipeline", "SeqArray"), \
                     vapply(c("GFisher", "GLOWr", "GLOWpipeline", "STAAR", "STAARpipeline", "SeqArray"), \
                            function(p) as.character(packageVersion(p)), ""))); \
  writeLines(lines, "/opt/glow-image-manifest.txt"); cat(lines, sep = "\n")'
