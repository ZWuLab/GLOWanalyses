#!/usr/bin/env Rscript
# ==============================================================================
# File: base-config/make-settings-reference.R
# Purpose: Generate base-config/settings-reference.md, the one table of every
#   GLOWanalyses setting for both routes, the `glow` command and the Terra workflow.
#   For each setting of the registry it gives the type, the default, the parts of the
#   workflow that read it, the Terra input that sets it, and the description. A
#   second table lists the inputs of cloud/terra/glow_scan.wdl with the setting each
#   one sets. The file is generated from three sources, so it cannot drift from them:
#     - base-config/settings-registry.yml       every setting, its type and default
#     - base-config/parts.yml                   which part reads which setting
#     - cloud/terra/workflow-settings.yml       which workflow input sets which setting
#   and from the input block and the parameter_meta block of the workflow.
#
# Usage (from the GLOWanalyses directory):
#   Rscript base-config/make-settings-reference.R            # write the file
#   Rscript base-config/make-settings-reference.R --check    # exit 1 if the file is stale
#
# The image build runs --check. With --check the last line on the screen is
# "Verdict for the settings reference: PASS" or "... FAIL".
suppressMessages(library(GLOWpipeline))

args <- commandArgs(trailingOnly = TRUE)
if (length(setdiff(args, "--check"))) stop("Usage: make-settings-reference.R [--check]")
check <- "--check" %in% args
OUT <- "base-config/settings-reference.md"
WDL <- "cloud/terra/glow_scan.wdl"
for (f in c("base-config/base_config.yml", "base-config/settings-registry.yml", "base-config/parts.yml",
            "cloud/terra/workflow-settings.yml", WDL))
  if (!file.exists(f)) stop("Run from the GLOWanalyses directory: ", f, " was not found.")

#################### THE SOURCES ####################

base <- read_glow_config("base-config/base_config.yml", base = "base-config/base_config.yml",
                         registry = "base-config/settings-registry.yml", root = ".")
entries <- attr(base, "glow_registry")$entries
parts <- glow_read_parts("base-config/parts.yml", "base-config/settings-registry.yml")
wf <- yaml::read_yaml("cloud/terra/workflow-settings.yml")

# The setting behind a name of the table of parts or of the workflow file: a field of a
# map (snv_set.calibration.report) belongs to the setting snv_set.calibration.
setting_of <- function(key) {
  k <- strsplit(key, ".", fixed = TRUE)[[1]]
  if (length(k) == 1L) key else paste(k[1:2], collapse = ".")
}
# The parts that name each setting, in any role: read, forwarded, an outside file, the
# path of an input, or a condition.
named_by <- function(p) {
  cond <- function(x) unlist(lapply(x, `[[`, "setting"))
  unique(vapply(c(p$settings, unlist(p$region_settings), p$forwarded, cond(p$run_when),
                  unlist(lapply(p$files, function(f) c(f$setting, f$pattern, cond(f$when)))),
                  unlist(lapply(p$inputs, function(i) c(i$setting, i$pattern, i$units_from, i$partial, cond(i$when))))),
                setting_of, character(1)))
}
readers <- list()
for (p in parts) for (k in named_by(p)) readers[[k]] <- c(readers[[k]], paste(p$step, p$part, sep = "/"))

# The workflow inputs that set each setting, and the settings a task writes on its own.
wf_inputs <- list()
for (n in names(wf$inputs)) for (k in c(wf$inputs[[n]]$setting, unlist(wf$inputs[[n]]$also))) {
  field <- sub(paste0("^", gsub(".", "\\.", setting_of(k), fixed = TRUE), "\\.?"), "", k)
  wf_inputs[[setting_of(k)]] <- c(wf_inputs[[setting_of(k)]], if (nzchar(field)) sprintf("`%s` (%s)", n, field) else sprintf("`%s`", n))
}
wf_task <- wf$task_values
names(wf_task) <- vapply(names(wf_task), setting_of, character(1))
# the sections whose settings the workflow's tasks read at all
WORKFLOW_SECTIONS <- c("top", "shared", "single_variant", "snv_set", "run")

# The inputs of the workflow, with their descriptions.
src <- readLines(WDL, warn = FALSE)
block <- function(open_pattern) {
  from <- grep("^workflow GLOWScan \\{", src)
  a <- from + grep(open_pattern, src[(from + 1L):length(src)])[1L]
  b <- a + grep("^  \\}", src[(a + 1L):length(src)])[1L]
  src[(a + 1L):(b - 1L)]
}
decl <- regmatches(block("^  input \\{"), regexec("^\\s+([A-Za-z\\[\\]]+\\??)\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*(=\\s*(.*))?$", block("^  input \\{"), perl = TRUE))
decl <- Filter(function(m) length(m) > 0L, decl)
meta <- regmatches(block("^  parameter_meta \\{"), regexec("^\\s+([A-Za-z_][A-Za-z0-9_]*):\\s*\"(.*)\"\\s*$", block("^  parameter_meta \\{")))
meta <- Filter(function(m) length(m) > 0L, meta)
meta_of <- stats::setNames(vapply(meta, `[[`, "", 3L), vapply(meta, `[[`, "", 2L))

#################### THE TEXT ####################

cell <- function(x) gsub("|", "\\|", gsub("\\s+", " ", x), fixed = TRUE)
code <- function(x) paste0("`", gsub("`", "'", x), "`")
type_words <- function(spec) {
  words <- switch(spec$type,
    logical = "true or false", integer = "whole number", number = "number", string = "text", path = "path",
    file_name = "file name", integer_list = "list of whole numbers", number_list = "list of numbers",
    string_list = "list of text", named_vector = "mapping of names to values", map = "mapping", sequence = "list",
    one_of = "one of several forms", r_function = "R function", spec$type)
  if (!is.null(spec$choices))
    words <- paste0(if (grepl("_list$", spec$type)) "list of: " else "one of ", paste(spec$choices, collapse = ", "))
  if (isTRUE(spec$null_ok)) words <- paste(words, "or null")
  words
}
default_words <- function(v) {
  if (is.null(v)) return("null")
  if (is.function(v)) return("an R function, see `base_config.yml`")
  if (is.list(v)) return(if (length(v)) "a mapping, see `base_config.yml`" else "empty")
  text <- if (is.logical(v)) tolower(as.character(v)) else as.character(v)
  one <- if (length(v) == 1L) text else paste0("[", paste(text, collapse = ", "), "]")
  if (nchar(one) > 70L) "a list, see `base_config.yml`" else code(one)
}
workflow_words <- function(e) {
  if (!is.null(wf_inputs[[e$key]])) return(paste(wf_inputs[[e$key]], collapse = ", "))
  if (!is.null(wf_task[[e$key]])) return(paste("written by the task:", wf_task[[e$key]]))
  if (e$section %in% WORKFLOW_SECTIONS) "the default" else "not run by the workflow"
}
SECTIONS <- c(top = "The three names of a run", data = "The raw inputs of a dataset (`data`)",
              shared = "Settings that several steps read (`shared`)", data_prep = "The data preparation (`data_prep`)",
              training = "The training of B and PI (`training`)", single_variant = "The single-variant step (`single_variant`)",
              snv_set = "The SNV-set step (`snv_set`)", run = "How the work runs (`run`)")
unknown <- setdiff(unique(vapply(entries, `[[`, "", "section")), names(SECTIONS))
if (length(unknown)) stop("The registry has a section this script does not know: ", paste(unknown, collapse = ", "))

lines <- c(
  "# GLOWanalyses settings reference", "",
  "Generated by `base-config/make-settings-reference.R` from the settings registry",
  "(`base-config/settings-registry.yml`), the table of parts (`base-config/parts.yml`), and the",
  "workflow with its declaration file (`cloud/terra/glow_scan.wdl`,",
  "`cloud/terra/workflow-settings.yml`). Do not edit it by hand. Regenerate it from the",
  "GLOWanalyses directory with `Rscript base-config/make-settings-reference.R`. The image build",
  "runs `--check` and fails when this file is stale.", "",
  "A configuration is a YAML file that sets some of these settings and takes the rest from",
  "`base-config/base_config.yml`. A setting is written under its section, for example",
  "`snv_set:` and below it `ld_threshold: 0.8`. The three names of a run have no section. The",
  "columns:", "",
  "- **Default**: the value of `base_config.yml`. A null means that the setting is off, or that",
  "  the workflow fills it in from the layout of the run (the description says which).",
  "- **Read by**: the parts whose completion record depends on the setting, as `step/part`",
  "  (`base-config/parts.yml`). `glow run <step>` runs the parts of a step in order.",
  "- **Terra input**: the input of `glow_scan.wdl` that sets the setting. \"the default\" means",
  "  that no input sets it, so the workflow runs the default. \"written by the task\" means that",
  "  a task of the workflow writes its own value, with the reason. \"not run by the workflow\"",
  "  marks the steps that the workflow does not run (the data preparation and the training).", "")
n_rows <- 0L
for (s in names(SECTIONS)) {
  of_section <- Filter(function(e) identical(e$section, s), entries)
  if (!length(of_section)) next
  lines <- c(lines, paste("##", SECTIONS[[s]]), "",
             "| Setting | Type | Default | Read by | Terra input | Description |", "|---|---|---|---|---|---|")
  for (e in of_section) {
    lines <- c(lines, sprintf("| %s | %s | %s | %s | %s | %s |", code(e$name), type_words(e$spec), default_words(e$default),
                              if (is.null(readers[[e$key]])) "" else paste(code(unique(readers[[e$key]])), collapse = ", "),
                              cell(workflow_words(e)), cell(e$description)))
    n_rows <- n_rows + 1L
  }
  lines <- c(lines, "")
}
lines <- c(lines, "## The inputs of the Terra workflow (`cloud/terra/glow_scan.wdl`)", "",
           "Every input with its type, its default in the workflow, the setting it sets, and its",
           "description. A default that differs from the default of the setting on purpose is marked",
           "with its reason. `cloud/terra/check-workflow-defaults.R` checks this table against the registry.", "",
           "| Input | Type | Default | Setting | Description |", "|---|---|---|---|---|")
for (m in decl) {
  n <- m[3L]; w <- wf$inputs[[n]]
  setting <- if (is.null(w)) "" else if (!is.null(w$no_setting)) paste("none:", w$no_setting)
             else paste0(paste(code(c(w$setting, unlist(w$also))), collapse = ", "),
                         if (!is.null(w$differs)) paste0(". The default differs: ", w$differs) else "")
  lines <- c(lines, sprintf("| %s | %s | %s | %s | %s |", code(n), code(m[2L]), if (nzchar(m[4L])) code(trimws(m[5L])) else "none",
                            cell(setting), cell(if (is.na(meta_of[n])) "" else meta_of[[n]])))
}
lines <- c(lines, "")

#################### WRITE OR CHECK ####################

if (check) {
  same <- file.exists(OUT) && identical(readLines(OUT, warn = FALSE), lines)
  if (!same) cat(OUT, "is not what the registry, the table of parts and the workflow give now.",
                 "Regenerate it: Rscript base-config/make-settings-reference.R\n")
  cat(sprintf("Verdict for the settings reference: %s\n", if (same) "PASS" else "FAIL"))
  quit(save = "no", status = if (same) 0L else 1L)
}
writeLines(lines, OUT)
cat(sprintf("Wrote %s: %d settings in %d sections, %d workflow inputs.\n", OUT, n_rows, length(unique(vapply(entries, `[[`, "", "section"))), length(decl)))
