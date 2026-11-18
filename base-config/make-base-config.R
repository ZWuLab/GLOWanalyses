#!/usr/bin/env Rscript
# ==============================================================================
# File: base-config/make-base-config.R
# Purpose: generate base-config/base_config.yml, the configuration with every
#   GLOWanalyses setting at its default, from base-config/settings-registry.yml.
#   The registry is the single source; this file is what a user reads and what
#   GLOWpipeline::read_glow_config() reads for the defaults. Run from the
#   GLOWanalyses root. Needs the GLOWpipeline package (the writer lives there).
#
# Usage:
#   Rscript base-config/make-base-config.R            # write the file
#   Rscript base-config/make-base-config.R --check    # exit 1 if the committed file is stale

args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--check" %in% args) "check" else "write"
REGISTRY <- "base-config/settings-registry.yml"
OUT_FILE <- "base-config/base_config.yml"
if (!file.exists(REGISTRY))
  stop("Run from the GLOWanalyses root (the directory holding base-config/).")

# The lines are a function of the registry alone (no date, no path of this machine),
# so two runs give the same bytes and a comparison with the committed file is exact.
generated <- GLOWpipeline::write_base_config(REGISTRY, file = NULL)

if (mode == "write") {
  writeLines(enc2utf8(generated), OUT_FILE, useBytes = TRUE)
  cat(sprintf("wrote %s: %d lines\n", OUT_FILE, length(generated)))
} else {
  if (!file.exists(OUT_FILE)) {
    cat("STALE: ", OUT_FILE, " does not exist; run the generator\n", sep = "")
    quit(status = 1L)
  }
  current <- readLines(OUT_FILE, warn = FALSE, encoding = "UTF-8")
  if (identical(current, generated)) {
    cat("base_config.yml is current\n")
  } else {
    n <- max(length(current), length(generated))
    pad <- function(x) c(x, rep("", n - length(x)))
    d <- which(pad(current) != pad(generated))
    cat(sprintf("STALE: %s differs from the generated text at %d line(s); first: line %d\n  committed: %s\n  generated: %s\n",
                OUT_FILE, length(d), d[1L],
                if (d[1L] <= length(current)) current[d[1L]] else "<missing>",
                if (d[1L] <= length(generated)) generated[d[1L]] else "<missing>"))
    quit(status = 1L)
  }
}
