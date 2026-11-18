#!/usr/bin/env Rscript
# ==============================================================================
# File: cloud/terra/check-workflow-defaults.R
# Purpose: Check that the Terra workflow and the settings registry agree. The two
#   routes of GLOWanalyses, `glow run` and cloud/terra/glow_scan.wdl, have to give
#   the same results for the same settings, so the workflow must not carry a default
#   of its own without saying so. The check compares
#     - the default of every workflow input with the default of the setting it sets;
#     - every value that a task writes into its configuration on its own with the
#       default of that setting;
#   through cloud/terra/workflow-settings.yml, which names the setting of each input
#   and declares each deliberate difference with its reason. A difference that is not
#   declared fails the check, and so does an input or a written setting that the file
#   does not know. The check also keeps the rule that input names stay distinct when
#   letter case is ignored (Terra's configuration store).
#
# Usage (from the GLOWanalyses directory):
#   Rscript cloud/terra/check-workflow-defaults.R [--wdl <file>] [--out <report.md>]
#
#   --wdl   the workflow to check. Default: cloud/terra/glow_scan.wdl.
#   --out   write the report as Markdown. Default: the screen.
#
# Output: a table per input and per task value, the counts, and as the last line
#   "Verdict for the workflow defaults: PASS" or "... FAIL". Exit status 1 on FAIL.
suppressMessages(library(GLOWpipeline))

args <- commandArgs(trailingOnly = TRUE)
opt <- list(wdl = "cloud/terra/glow_scan.wdl", out = NA_character_)
i <- 1L
while (i <= length(args)) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opt) || i == length(args)) stop("Usage: check-workflow-defaults.R [--wdl <file>] [--out <report.md>]")
  opt[[key]] <- args[i + 1L]; i <- i + 2L
}
for (f in c("base-config/base_config.yml", "base-config/settings-registry.yml", "cloud/terra/workflow-settings.yml", opt$wdl))
  if (!file.exists(f)) stop("Run from the GLOWanalyses directory: ", f, " was not found.")

#################### THE THREE SOURCES ####################

# The defaults of the registry, typed by the loader.
base <- read_glow_config("base-config/base_config.yml", base = "base-config/base_config.yml",
                         registry = "base-config/settings-registry.yml", root = ".")
registered <- vapply(attr(base, "glow_registry")$entries, function(e)
  if (e$section == "top") e$name else paste(e$section, e$name, sep = "."), character(1))
default_of <- function(key) { v <- base; for (n in strsplit(key, ".", fixed = TRUE)[[1]]) { if (!is.list(v)) return(NULL); v <- v[[n]] }; v }
is_registered <- function(key) paste(utils::head(strsplit(key, ".", fixed = TRUE)[[1]], 2L), collapse = ".") %in% registered || key %in% registered

map <- yaml::read_yaml("cloud/terra/workflow-settings.yml")

# The inputs of the workflow: type, name and default, from the input block of GLOWScan.
src <- readLines(opt$wdl, warn = FALSE)
wf <- grep("^workflow GLOWScan \\{", src)
open <- wf + grep("^  input \\{", src[(wf + 1L):length(src)])[1L]
close <- open + grep("^  \\}", src[(open + 1L):length(src)])[1L]
decl <- regmatches(src[(open + 1L):(close - 1L)],
                   regexec("^\\s+([A-Za-z\\[\\]]+\\??)\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*(=\\s*(.*))?$", src[(open + 1L):(close - 1L)], perl = TRUE))
decl <- Filter(function(m) length(m) > 0L, decl)
inputs <- data.frame(name = vapply(decl, `[[`, "", 3L), type = vapply(decl, `[[`, "", 2L),
                     default = vapply(decl, function(m) if (nzchar(m[4L])) trimws(m[5L]) else NA_character_, ""),
                     stringsAsFactors = FALSE)
# A WDL literal as an R value: a number, true or false, a string, or a list of strings.
wdl_value <- function(text, type) {
  if (grepl("^Array", type)) return(unlist(regmatches(text, gregexpr("\"[^\"]*\"", text))) |> gsub(pattern = "^\"|\"$", replacement = ""))
  if (type == "Boolean") return(identical(text, "true"))
  if (type == "Int") return(as.integer(text))
  if (type == "Float") return(as.numeric(text))
  gsub("^\"|\"$", "", text)
}
same <- function(a, b) {
  if (is.null(a) || is.null(b)) return(is.null(a) && is.null(b))
  if (is.numeric(a) && is.numeric(b)) return(length(a) == length(b) && all(a == b))
  identical(unname(a), unname(b))
}
show <- function(v) if (is.null(v)) "null" else substr(paste(format(unlist(v)), collapse = ", "), 1L, 60L)

#################### THE INPUTS ####################

rows <- lapply(seq_len(nrow(inputs)), function(k) {
  n <- inputs$name[k]; m <- map$inputs[[n]]
  row <- list(input = n, setting = "", workflow = if (is.na(inputs$default[k])) "(none)" else inputs$default[k], registry = "", verdict = "")
  if (is.null(m)) { row$verdict <- "FAIL: not in workflow-settings.yml"; return(row) }
  if (!is.null(m$no_setting)) { row$verdict <- paste("no setting:", m$no_setting); return(row) }
  keys <- c(m$setting, unlist(m$also))
  row$setting <- paste(keys, collapse = ", ")
  unknown <- keys[!vapply(keys, is_registered, logical(1))]
  if (length(unknown)) { row$verdict <- paste("FAIL: not a registered setting:", paste(unknown, collapse = ", ")); return(row) }
  row$registry <- paste(vapply(keys, function(key) show(default_of(key)), character(1)), collapse = " | ")
  if (is.na(inputs$default[k])) { row$verdict <- "no default in the workflow"; return(row) }
  value <- wdl_value(inputs$default[k], inputs$type[k])
  if (!is.null(m$values)) {
    if (!value %in% names(m$values)) { row$verdict <- "FAIL: the default is not among the declared words"; return(row) }
    value <- m$values[[value]]
  }
  equal <- all(vapply(keys, function(key) same(value, default_of(key)), logical(1)))
  row$verdict <- if (equal && is.null(m$differs)) "same default"
                 else if (!equal && !is.null(m$differs)) paste("differs, declared:", m$differs)
                 else if (equal) "FAIL: declared as differing, and the defaults are the same"
                 else "FAIL: the defaults differ, and no reason is declared"
  row
})
tab <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
stale <- setdiff(names(map$inputs), inputs$name)
lower <- tolower(inputs$name)
clash <- inputs$name[lower %in% lower[duplicated(lower)]]

#################### THE VALUES THAT TASKS WRITE ####################

# Every configuration that a task writes is a block between `<<'YCONF'` and `YCONF`,
# plus the lines that a task appends with echo or printf. A line `key: value` is read
# with its section. A line whose value holds a placeholder takes the value of an
# input, so its key has to be the setting of an input. A line with a literal value is
# a value the task writes on its own: equal to the registry default, or declared.
starts <- grep("<<'YCONF'", src, fixed = TRUE); ends <- grep("^\\s*YCONF\\s*$", src)
input_keys <- unique(unlist(lapply(map$inputs, function(m) c(m$setting, unlist(m$also)))))
# a map that a task builds from inputs (the variant filter) counts through its fields
from_input <- function(key) key %in% input_keys || any(startsWith(input_keys, paste0(key, ".")))
task_rows <- list(); map_rows <- list()
# The maps that a block writes on one line, `key: {field: ..., field: ...}`, wherever the line
# stands: also inside the text of a conditional placeholder. Each field has to be a field of
# the setting's map in the registry.
inline_maps <- function(body, section, line) {
  found <- regmatches(body, gregexpr("[A-Za-z_]+: \\{[^{}]*(~\\{[^}]*\\}[^{}]*)*\\}", body))[[1L]]
  lapply(found, function(text) {
    key <- sub(":.*$", "", text)
    inner <- sub("^[A-Za-z_]+: \\{", "", sub("\\}$", "", text))
    # the fields: split at each comma that a field name and a colon follow (a name may be quoted)
    pairs <- strsplit(inner, ",(?=\\s*[\"']?[A-Za-z_]+[\"']?\\s*:)", perl = TRUE)[[1L]]
    fields <- lapply(pairs, function(p) {
      m <- regmatches(p, regexec("^\\s*[\"']?([A-Za-z_]+)[\"']?\\s*:\\s*(.*)$", p))[[1L]]
      if (length(m)) list(field = m[2L], value = trimws(m[3L])) else list(field = NA_character_, value = trimws(p))
    })
    list(key = if (nzchar(section)) paste(section, key, sep = ".") else key, fields = fields, line = line)
  })
}
for (b in seq_along(starts)) {
  block <- src[(starts[b] + 1L):(ends[ends > starts[b]][1L] - 1L)]
  indent <- min(nchar(sub("\\S.*$", "", block[nzchar(trimws(block))])))
  section <- ""
  for (line in block) {
    body <- substring(line, indent + 1L)
    if (!nzchar(trimws(body)) || grepl("^\\s*#", body)) next
    map_rows <- c(map_rows, inline_maps(body, section, starts[b]))
    m <- regmatches(body, regexec("^( *)([A-Za-z_]+):\\s*(.*)$", body))[[1L]]
    if (length(m) == 0L) next                                  # the body of a block scalar, or a placeholder line
    depth <- nchar(m[2L]); key <- m[3L]; value <- sub("\\s+#.*$", "", m[4L])
    if (depth == 0L && !nzchar(value)) { section <- key; next }
    if (depth > 2L) next
    full <- if (depth == 0L) key else paste(section, key, sep = ".")
    task_rows[[length(task_rows) + 1L]] <- list(key = full, value = value, line = starts[b])
  }
}
# the lines appended after a block: `  n_cores: ${N_CORES}` under the last section, which is `shared`
appended <- grep("^\\s*(echo|printf) .*n_cores: ", src)
for (a in appended) task_rows[[length(task_rows) + 1L]] <- list(key = "shared.n_cores", value = "${N_CORES}", line = a)
# The B function is not written in a block: the task gives the text of the input `b_func`, in a
# file that the engine writes (write_lines), to cloud/terra/write-b-func.R, which writes the
# setting. The call counts as the setting taking the input when the file is that text.
text_decl <- any(grepl("^\\s*File\\s+b_func_text\\s*=\\s*write_lines\\(\\[b_func\\]\\)\\s*$", src))
for (a in grep("^[^#]*write-b-func\\.R", src))
  task_rows[[length(task_rows) + 1L]] <- list(key = "snv_set.b_func", line = a,
    value = if (text_decl && grepl("'~{b_func_text}'", src[a], fixed = TRUE)) "~{b_func}" else "\"a file that is not the text of the input b_func\"")
# The inputs of the workflow that set a key, or a field of it.
inputs_of_key <- function(key) names(Filter(function(m) {
  keys <- as.character(c(m$setting, unlist(m$also))); length(keys) > 0L && any(keys == key | startsWith(keys, paste0(key, ".")))
}, map$inputs))
# The names that the placeholders of a value use, apart from those inside defined(...) and
# inside quoted text: the condition of a branch may name another input, its value may not.
names_used <- function(value) {
  inner <- unlist(regmatches(value, gregexpr("~\\{[^}]*\\}", value)))
  inner <- gsub("\"[^\"]*\"|'[^']*'", "", gsub("defined\\([^)]*\\)", "", inner))
  unique(unlist(regmatches(inner, gregexpr("[A-Za-z_][A-Za-z0-9_]*", inner))))
}
# What the loader says to one setting with a literal value: "" when it reads it.
loader_says <- function(key, text) {
  k <- strsplit(key, ".", fixed = TRUE)[[1L]]
  doc <- if (length(k) == 1L) sprintf("%s: %s", key, text) else c(paste0(k[1L], ":"), sprintf("  %s: %s", k[2L], text))
  f <- tempfile(fileext = ".yml"); on.exit(unlink(f, expand = FALSE)); writeLines(doc, f)
  tryCatch({ read_glow_config(f, base = "base-config/base_config.yml", registry = "base-config/settings-registry.yml", root = "."); "" },
           error = function(e) trimws(gsub("\\s+", " ", gsub(f, "", sub("^.*problem\\(s\\):", "", conditionMessage(e)), fixed = TRUE))))
}
task_tab <- unique(do.call(rbind, lapply(task_rows, function(r) {
  placeholder <- grepl("~{", r$value, fixed = TRUE) || grepl("${", r$value, fixed = TRUE)
  declared <- map$task_values[[r$key]]
  # an input of the workflow in the value of a setting that it does not set is a crossed wire
  stray <- if (placeholder && is_registered(r$key)) setdiff(intersect(names_used(r$value), inputs$name), inputs_of_key(r$key)) else character(0)
  # a literal value has to be one the loader reads for this setting, declared or not
  refused <- if (!placeholder && is_registered(r$key) && length(strsplit(r$key, ".", fixed = TRUE)[[1L]]) <= 2L) loader_says(r$key, r$value) else ""
  verdict <- if (!is_registered(r$key)) "FAIL: not a registered setting"
    else if (length(stray)) paste0("FAIL: takes the input ", paste(stray, collapse = ", "), ", which workflow-settings.yml does not give to this setting")
    else if (nzchar(refused)) paste0("FAIL: the loader refuses the value: ", substr(refused, 1L, 160L))
    else if (placeholder && from_input(r$key)) "the value of an input"
    else if (!is.null(declared)) paste("written by the task, declared:", declared)
    else if (placeholder) "FAIL: takes a value that no input of workflow-settings.yml sets"
    else if (same(yaml::yaml.load(r$value), default_of(r$key))) "written by the task, same as the registry default"
    else "FAIL: written by the task, differs from the registry default, and no reason is declared"
  data.frame(setting = r$key, written = substr(r$value, 1L, 50L), registry = show(default_of(r$key)), verdict = verdict, stringsAsFactors = FALSE)
})))
never_written <- setdiff(names(map$task_values), task_tab$setting)
# The fields of the maps written on one line, each on its own. A field has to be a field of the
# setting. A field whose value takes an input (a placeholder, or a name joined into the text of a
# conditional line) may take only an input that workflow-settings.yml gives to that field: two
# inputs written into each other's fields are a crossed wire. A field with a literal value is
# given to the loader, alone in its map, and fails when the loader refuses it.
map_tab <- unique(do.call(rbind, lapply(map_rows, function(r) {
  known <- names(default_of(r$key))
  problems <- character(0)
  for (f in r$fields) {
    if (is.na(f$field)) { problems <- c(problems, paste0("not a `field: value` pair: ", substr(f$value, 1L, 40L))); next }
    if (is.null(known) || !f$field %in% known) { problems <- c(problems, paste0("not a field of the setting: ", f$field)); next }
    joined <- grepl("['\"]\\s*\\+|\\+\\s*['\"(]", f$value)                # text joined with a name in a conditional line
    takes_input <- grepl("~{", f$value, fixed = TRUE) || joined
    if (takes_input) {
      used <- if (joined) unique(unlist(regmatches(gsub("defined\\([^)]*\\)", "", f$value), gregexpr("[A-Za-z_][A-Za-z0-9_]*", gsub("defined\\([^)]*\\)", "", f$value)))))
              else names_used(f$value)
      stray <- setdiff(intersect(used, inputs$name), inputs_of_key(paste(r$key, f$field, sep = ".")))
      if (length(stray)) problems <- c(problems, paste0("the field ", f$field, " takes the input ", paste(stray, collapse = ", "),
                                                        ", which workflow-settings.yml does not give to it"))
    } else if (length(strsplit(r$key, ".", fixed = TRUE)[[1L]]) <= 2L) {
      refused <- loader_says(r$key, sprintf("{%s: %s}", f$field, f$value))
      if (nzchar(refused)) problems <- c(problems, paste0("the loader refuses the field ", f$field, ": ", substr(refused, 1L, 120L)))
    }
  }
  data.frame(setting = r$key, fields = paste(vapply(r$fields, function(f) if (is.na(f$field)) "?" else f$field, character(1)), collapse = ", "),
             verdict = if (!is_registered(r$key)) "FAIL: not a registered setting"
                       else if (is.null(known)) "FAIL: the setting is not a map with fields in the registry"
                       else if (length(problems)) paste0("FAIL: ", paste(problems, collapse = "; "))
                       else "every field is a field of the setting, with its own input or a value the loader reads", stringsAsFactors = FALSE)
})))
if (is.null(map_tab)) map_tab <- data.frame(setting = character(0), fields = character(0), verdict = character(0))

#################### THE REPORT ####################

kind <- function(v) sub(":.*$", "", sub(",.*$", "", v))
n_fail <- sum(grepl("^FAIL", tab$verdict)) + sum(grepl("^FAIL", task_tab$verdict)) + sum(grepl("^FAIL", map_tab$verdict)) +
  length(stale) + length(clash) + length(never_written)
pass <- n_fail == 0L && nrow(tab) > 0L && nrow(task_tab) > 0L
md_table <- function(d) c(paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("| ", paste(rep("---", ncol(d)), collapse = " | "), " |"),
                          apply(d, 1L, function(r) paste0("| ", paste(gsub("\\|", "/", r), collapse = " | "), " |")))
with_default <- !is.na(inputs$default)
lines <- c(
  "# The Terra workflow against the settings registry", "",
  sprintf("Written by `cloud/terra/check-workflow-defaults.R` for `%s`.", opt$wdl), "",
  "## Summary", "",
  sprintf("- Inputs of the workflow: %d, of which %d declare a default.", nrow(inputs), sum(with_default)),
  sprintf("- Defaults that are the default of their setting: %d.", sum(tab$verdict == "same default")),
  sprintf("- Defaults that differ, each declared with its reason: %d.", sum(grepl("^differs, declared", tab$verdict))),
  sprintf("- Inputs without a setting, each with its reason: %d.", sum(grepl("^no setting", tab$verdict))),
  sprintf("- Inputs without a default in the workflow: %d.", sum(tab$verdict == "no default in the workflow")),
  sprintf("- Settings that a task writes: %d, of which %d are the value of an input, %d are declared values of the task and %d equal the registry default.",
          nrow(task_tab), sum(task_tab$verdict == "the value of an input"), sum(grepl("^written by the task, declared", task_tab$verdict)),
          sum(task_tab$verdict == "written by the task, same as the registry default")),
  sprintf("- Failing: %d.", n_fail),
  sprintf("- Verdict: **%s**.", if (pass) "PASS" else "FAIL"), "",
  if (length(stale)) c(sprintf("Entries of workflow-settings.yml for inputs that the workflow does not have: %s.", paste(stale, collapse = ", ")), ""),
  if (length(clash)) c(sprintf("Input names that are the same when letter case is ignored: %s.", paste(clash, collapse = ", ")), ""),
  if (length(never_written)) c(sprintf("Task values declared in workflow-settings.yml that no task writes: %s.", paste(never_written, collapse = ", ")), ""),
  "## Inputs", "", md_table(tab), "",
  "## Settings that the tasks write", "", md_table(task_tab), "",
  "## Maps that the tasks write on one line", "", if (nrow(map_tab)) md_table(map_tab) else "None.", "")
if (!is.na(opt$out)) { writeLines(lines, opt$out); cat("Report written to", opt$out, "\n") } else cat(lines, sep = "\n")
if (n_fail > 0L) {
  cat("\nFailing:\n")
  print(rbind(tab[grepl("^FAIL", tab$verdict), c("input", "verdict")] |> stats::setNames(c("what", "verdict")),
              task_tab[grepl("^FAIL", task_tab$verdict), c("setting", "verdict")] |> stats::setNames(c("what", "verdict")),
              map_tab[grepl("^FAIL", map_tab$verdict), c("setting", "verdict")] |> stats::setNames(c("what", "verdict"))), row.names = FALSE)
}
cat(sprintf("Inputs: %d. Task settings: %d. Maps on one line: %d. Failing: %d.\n", nrow(inputs), nrow(task_tab), nrow(map_tab), n_fail))
cat(sprintf("Verdict for the workflow defaults: %s\n", if (pass) "PASS" else "FAIL"))
quit(save = "no", status = if (pass) 0L else 1L)
