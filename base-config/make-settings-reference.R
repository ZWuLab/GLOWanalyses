#!/usr/bin/env Rscript
# ==============================================================================
# File: base-config/make-settings-reference.R
# Purpose: generate base-config/settings-reference.md, the one table of every GLOWanalyses
#   setting for both routes: each top-level field of the four base configuration files with
#   its default, the layer a run sets it in (the shipped example's _cohort.R and run configs),
#   the templates that read it and the Terra workflow input that carries it; the fields the
#   templates read with an in-code default; and every input of cloud/terra/glow_scan.wdl with
#   its type, default and description. Base R only, so the image and the release gate can
#   run it. Run from the GLOWanalyses root.
#
# Usage:
#   Rscript base-config/make-settings-reference.R              # write the file
#   Rscript base-config/make-settings-reference.R --check      # exit 1 if the committed file is stale
#   Rscript base-config/make-settings-reference.R --self-test  # exit 1 if an invariant fails
#
#   - a base setting is a top-level `name <- value` of a base file; its description is the
#     comment on the assignment's last line (plus indented continuation lines), else the
#     contiguous comment block right above it, else a `# ---- title ----` header line;
#   - a required setting is one the base files document in a `#   name  description` list
#     without assigning it, which the example's _cohort.R or run configs set;
#   - "set by" marks the example's layers (cohort, run), else base;
#   - "read by" is a word match on the templates' code, comment lines removed;
#   - the workflow column follows a `~{input}` placeholder from a task's config.R block to the
#     task input and through the workflow's `call` block to the workflow input.

args <- commandArgs(trailingOnly = TRUE)
mode <- if ("--check" %in% args) "check" else if ("--self-test" %in% args) "self-test" else "write"
if (!file.exists("base-config/snv_set_base.R") || !file.exists("cloud/terra/glow_scan.wdl"))
  stop("Run from the GLOWanalyses root (the directory holding base-config/ and cloud/).")

OUT_FILE <- "base-config/settings-reference.md"
WDL_FILE <- "cloud/terra/glow_scan.wdl"

# The four stages: base file, title, the template directories whose scripts read the stage's
# settings (04-summary reads the SNV-set configuration), and the example run configurations.
STAGES <- list(
  data_prep      = list(file = "base-config/data_prep_base.R",      title = "Stage 00, data preparation",
                        dirs = "00-data-prep",                        runs = "runs/example/data-prep/config.R"),
  training       = list(file = "base-config/training_base.R",       title = "Stage 01, training",
                        dirs = "01-training",                         runs = "runs/example/training/config.R"),
  single_variant = list(file = "base-config/single_variant_base.R", title = "Stage 02, single-variant scan and calibration",
                        dirs = "02-single-variant",                   runs = "runs/example/single-variant/config.R"),
  snv_set        = list(file = "base-config/snv_set_base.R",        title = "Stages 03 and 04, SNV-set scan and summary",
                        dirs = c("03-snv-set", "04-summary"),
                        runs = c("runs/example/snv-set/config.R", "runs/example/snv-set-window/config.R",
                                 "runs/example/snv-set-coding/config.R")))
COHORT_FILE <- "runs/example/_cohort.R"

# Task-scoped names a config.R placeholder may use that are not task inputs, and the workflow
# input each stands for.
TASK_VARS <- c(chrom = "chromosomes", use_spa_r = "use_spa", N_CORES = "n_cores")
# Outputs of one task that another task's call passes on, and the workflow input behind each.
DERIVED <- c(chromosomes_checked = "chromosomes")

#################### THE BASE FILES ####################

# Does a comment line document a field other than `name` (a `# other : text` or `#   other  text`
# entry, or a commented-out `# other <- value` example)?
other_fields_doc <- function(line, name) {
  m <- regmatches(line, regexec("^#\\s*([A-Za-z_][A-Za-z0-9_.]*)\\s*(:|<-|\\s{2,}\\S)", line))[[1L]]
  length(m) > 0L && m[2L] != name && !m[2L] %in% c("e.g", "i.e", "NOTE", "Note")
}

# Strip the comment marks and join a block of comment lines into one description.
clean_comment <- function(lines) {
  if (!length(lines)) return("")
  x <- sub("^\\s*#+\\s?", "", lines)
  x <- gsub("\\s+", " ", trimws(paste(x, collapse = " ")))
  gsub("|", "\\|", x, fixed = TRUE)
}

# One row per top-level assignment of an R file: name, default (deparsed), description, line.
parse_assignments <- function(path) {
  src <- readLines(path, warn = FALSE)
  ex  <- parse(path, keep.source = TRUE)
  pd  <- utils::getParseData(ex, includeText = TRUE)
  src_line <- vapply(attr(ex, "srcref"), function(s) as.integer(s[1L]), integer(1))
  top <- pd[pd$parent == 0 & pd$token == "expr", ]
  top <- top[order(top$line1), ]
  comments <- pd[pd$token == "COMMENT", ]
  rows <- list()
  for (k in seq_len(nrow(top))) {
    e <- top[k, ]
    kids <- pd[pd$parent == e$id, ]
    kids <- kids[order(kids$line1, kids$col1), ]
    if (nrow(kids) < 3L || kids$token[2L] != "LEFT_ASSIGN") next
    lhs <- pd[pd$parent == kids$id[1L], ]
    if (!any(lhs$token == "SYMBOL")) next
    name <- lhs$text[lhs$token == "SYMBOL"][1L]
    call <- ex[[which(src_line == e$line1)[1L]]]
    default <- paste(deparse(call[[3L]], width.cutoff = 500L), collapse = " ")
    default <- gsub("\\s+", " ", default)
    # the inline comment on the last line, with indented continuation lines
    inl <- comments[comments$line1 == e$line2 & comments$col1 > 1L, ]
    if (nrow(inl)) {
      desc <- inl$text[1L]
      ln <- e$line2 + 1L
      while (ln <= length(src) && grepl("^\\s+#", src[ln])) { desc <- c(desc, src[ln]); ln <- ln + 1L }
    } else {
      # the contiguous comment block right above, stopping at a blank line or a rule line
      ln <- e$line1 - 1L; desc <- character()
      while (ln >= 1L && grepl("^#", src[ln]) && !grepl("^#\\s*[-=]{3,}", src[ln])) { desc <- c(src[ln], desc); ln <- ln - 1L }
      # a block may hold documentation entries (`# field : text`, `#   field  text`) for several
      # fields: keep this field's own entry with its continuation lines; a block that opens with
      # another field's entry or a commented-out example assignment is not this field's description
      if (length(desc)) {
        starts <- which(vapply(desc, function(l) length(regmatches(l, regexec("^#\\s*([A-Za-z_][A-Za-z0-9_.]*)\\s*(:|\\s{2,}\\S)", l))[[1L]]) > 0L, logical(1)))
        own <- starts[vapply(desc[starts], function(l) identical(sub("^#\\s*([A-Za-z_][A-Za-z0-9_.]*).*$", "\\1", l), name), logical(1))]
        if (length(own)) {
          nxt <- starts[starts > own[1L]]
          desc <- desc[own[1L]:(if (length(nxt)) nxt[1L] - 1L else length(desc))]
        } else if (other_fields_doc(desc[1L], name)) { ln <- ln + length(desc); desc <- character() }
      }
      if (!length(desc) && ln >= 1L && grepl("^#\\s*-{3,}\\s*\\S", src[ln])) {
        desc <- sub("^#\\s*-{3,}\\s*", "", src[ln]); desc <- sub("\\s*-{3,}\\s*$", "", desc)   # a `# ---- title ----` header
      }
      if (!length(desc)) {
        # the group's comment: walk up past the assignments of the same group, past comment
        # lines that document other fields and past indented continuation lines, to the nearest
        # top-level comment or header
        while (ln >= 1L && (!grepl("^#", src[ln]) || grepl("^#\\s{2,}", src[ln]) || other_fields_doc(src[ln], name))) ln <- ln - 1L
        if (ln >= 1L && grepl("^#\\s*-{3,}\\s*\\S", src[ln])) {
          desc <- sub("^#\\s*-{3,}\\s*", "", src[ln]); desc <- sub("\\s*-{3,}\\s*$", "", desc)
        } else {
          while (ln >= 1L && grepl("^#", src[ln]) && !grepl("^#\\s*[-=]{3,}", src[ln])) { desc <- c(src[ln], desc); ln <- ln - 1L }
        }
        if (length(desc)) desc <- c("(group)", desc)
      }
    }
    description <- clean_comment(desc)
    description <- sub(paste0("^", name, "\\s*:\\s*"), "", description)      # the entry's own name
    rows[[length(rows) + 1L]] <- data.frame(name = name, default = default,
                                            description = description, line = e$line1,
                                            stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(data.frame(name = character(), default = character(), description = character(), line = integer()))
  do.call(rbind, rows)
}

# The `#   name  description` documentation lists of a base file (the fields a run must set),
# with continuation lines indented past the name.
doc_list <- function(path) {
  src <- readLines(path, warn = FALSE)
  out <- list(); cur <- NULL
  for (l in src) {
    m <- regmatches(l, regexec("^#\\s{2,}([A-Za-z_][A-Za-z0-9_.]*)\\s{2,}(\\S.*)$", l))[[1L]]          # `#   name  text`
    if (!length(m)) m <- regmatches(l, regexec("^#\\s{2,}([A-Za-z]+_[A-Za-z0-9_.]*)\\s(\\S.*)$", l))[[1L]]   # one space, underscored name
    if (!length(m)) m <- regmatches(l, regexec("^#\\s*([A-Za-z_][A-Za-z0-9_.]*)\\s*:\\s*(\\S.*)$", l))[[1L]]   # `# name : text`
    if (!length(m)) {                                                                                   # `# name <- example   # text`
      mm <- regmatches(l, regexec("^#\\s*([A-Za-z_][A-Za-z0-9_.]*)\\s*<-\\s*[^#]*#\\s*(\\S.*)$", l))[[1L]]
      if (length(mm)) m <- mm
    }
    if (length(m) && !m[2L] %in% c("e.g", "i.e")) { cur <- m[2L]; if (is.null(out[[cur]])) out[[cur]] <- m[3L] else cur <- NULL; next }
    if (!is.null(cur) && grepl("^#\\s{2,}\\S", l) && !grepl("<-", l)) { out[[cur]] <- paste(out[[cur]], sub("^#\\s+", "", l)); next }
    cur <- NULL
  }
  vapply(out, function(x) clean_comment(paste0("# ", x)), character(1))
}

#################### THE TEMPLATES ####################

template_scripts <- function(dirs) {
  f <- list.files(dirs, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  f[!grepl("/_", f)]                        # a leading underscore marks a shared library file
}
code_of <- function(path) {
  x <- readLines(path, warn = FALSE)
  x[!grepl("^\\s*#", x)]
}
readers_of <- function(name, scripts, codes) {
  hit <- vapply(codes, function(cd) any(grepl(paste0("\\b", name, "\\b"), cd)), logical(1))
  scripts[hit]
}
fmt_readers <- function(r) {
  if (!length(r)) return("")
  r <- sub("^0[0-9]-[a-z-]+/", "", r)       # the stage directory is the table's heading
  if (length(r) > 4L) paste(c(sprintf("`%s`", r[1:4]), sprintf("and %d more", length(r) - 4L)), collapse = ", ")
  else paste(sprintf("`%s`", r), collapse = ", ")
}

# Optional settings the templates read through g0()/get0(): name, scripts, default expression.
incode_fields <- function(scripts, codes) {
  rows <- list()
  for (i in seq_along(scripts)) {
    for (l in codes[[i]]) {
      for (m in regmatches(l, gregexpr('(g0|get0)\\("([A-Za-z_][A-Za-z0-9_.]*)"', l))[[1L]]) {
        nm <- sub('^(g0|get0)\\("([^"]*)".*$', "\\2", m)
        dflt <- regmatches(l, regexec(paste0('null_or\\(\\s*g0\\("', nm, '"\\)\\s*,\\s*(.+)\\)\\s*(#.*)?$'), l))[[1L]]
        if (!length(dflt)) dflt <- regmatches(l, regexec(paste0('g0\\("', nm, '"\\s*,\\s*(.+)\\)\\s*(#.*)?$'), l))[[1L]]
        if (!length(dflt)) dflt <- regmatches(l, regexec(paste0('get0\\("', nm, '"[^)]*ifnotfound\\s*=\\s*([^,)]+)'), l))[[1L]]
        rows[[length(rows) + 1L]] <- data.frame(name = nm, script = scripts[i],
                                                default = if (length(dflt)) trimws(dflt[2L]) else "see the script",
                                                stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(data.frame(name = character(), script = character(), default = character()))
  r <- do.call(rbind, rows); r[!duplicated(r[c("name", "script")]), ]
}

#################### THE WORKFLOW ####################

# A block of the WDL from a line matching `open` to the first following line matching `close`.
block_lines <- function(src, open, close, from = 1L) {
  s <- grep(open, src); s <- s[s >= from][1L]
  if (is.na(s)) return(integer())
  e <- grep(close, src); e <- e[e > s][1L]
  (s + 1L):(e - 1L)
}

parse_wdl <- function(path) {
  src <- readLines(path, warn = FALSE)
  wf  <- grep("^workflow GLOWScan \\{", src)
  # workflow inputs: `Type name = default` or `Type name`
  inp <- src[block_lines(src, "^  input \\{", "^  \\}$", from = wf)]
  inputs <- do.call(rbind, lapply(inp[grepl("^\\s+[A-Za-z]", inp)], function(l) {
    m <- regmatches(l, regexec("^\\s+(\\S+)\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*(=\\s*(.*?))?\\s*$", l))[[1L]]
    data.frame(name = m[3L], type = m[2L], default = if (nzchar(m[4L])) m[5L] else "", stringsAsFactors = FALSE)
  }))
  # descriptions
  pm <- src[block_lines(src, "^  parameter_meta \\{", "^  \\}$", from = wf)]
  meta <- list()
  for (l in pm) { m <- regmatches(l, regexec('^\\s+([A-Za-z_][A-Za-z0-9_]*):\\s*"(.*)"\\s*$', l))[[1L]]; if (length(m)) meta[[m[2L]]] <- m[3L] }
  inputs$description <- vapply(inputs$name, function(n) if (is.null(meta[[n]])) "" else gsub("|", "\\|", meta[[n]], fixed = TRUE), character(1))
  # the call blocks: task input -> workflow expression
  calls <- list()
  for (s in grep("^\\s+call ([A-Za-z]+) \\{", src)) {
    task <- sub("^\\s+call ([A-Za-z]+) \\{.*$", "\\1", src[s])
    e <- s; depth <- 0L
    repeat { depth <- depth + lengths(regmatches(src[e], gregexpr("\\{", src[e]))) - lengths(regmatches(src[e], gregexpr("\\}", src[e]))); if (depth <= 0L) break; e <- e + 1L }
    txt <- paste(src[s:e], collapse = " ")
    txt <- sub("^.*input:\\s*", "", txt); txt <- sub("\\s*\\}\\s*$", "", txt)
    pos <- gregexpr("\\b[A-Za-z_][A-Za-z0-9_]*\\s*=\\s*", txt)[[1L]]
    names_ <- sub("\\s*=\\s*$", "", regmatches(txt, list(pos))[[1L]])
    ends <- c(pos[-1L] - 1L, nchar(txt))
    vals <- trimws(sub(",\\s*$", "", substring(txt, pos + attr(pos, "match.length"), ends)))
    calls[[task]] <- setNames(vals, names_)
  }
  # the tasks' config.R lines: field <- expression with placeholders
  task_starts <- grep("^task ([A-Za-z]+) \\{", src)
  task_names  <- sub("^task ([A-Za-z]+) \\{.*$", "\\1", src[task_starts])
  task_ends   <- c(task_starts[-1L] - 1L, length(src))
  maps <- list()
  for (t in seq_along(task_starts)) {
    body <- src[task_starts[t]:task_ends[t]]
    hd <- grep("cat > /tmp/glow/run/config.R <<'RCONF'", body)
    lines <- character()
    if (length(hd)) { en <- grep("^\\s*RCONF\\s*$", body); en <- en[en > hd[1L]][1L]; lines <- body[(hd[1L] + 1L):(en - 1L)] }
    lines <- c(lines, sub('^\\s*echo "([^"]*)".*$', "\\1", body[grepl('^\\s*echo "[A-Za-z_]+\\s*<-.*>> /tmp/glow/run/config.R', body)]))
    last <- NULL
    for (l in lines) {
      if (grepl("^\\s*(#|source\\()", l)) next
      fields <- sub("\\s*<-$", "", regmatches(l, gregexpr("\\b[A-Za-z_][A-Za-z0-9_.]*\\s*<-", l))[[1L]])
      if (!length(fields)) { if (is.null(last)) next; fields <- last }
      last <- fields[length(fields)]
      ph <- unlist(regmatches(l, gregexpr("~\\{[^}]*\\}|\\$\\{[A-Z_]+\\}", l)))
      ids <- unique(unlist(lapply(ph, function(p) {
        p <- gsub("'[^']*'|\"[^\"]*\"", " ", p)          # drop the quoted literals of a conditional
        w <- regmatches(p, gregexpr("[A-Za-z_][A-Za-z0-9_]*", p))[[1L]]
        setdiff(w, c("if", "then", "else", "defined", "select_first", "sep", "true", "false"))
      })))
      for (f in fields) {
        res <- character()
        for (id in ids) {
          if (!is.null(TASK_VARS[id]) && !is.na(TASK_VARS[id])) { res <- c(res, TASK_VARS[[id]]); next }
          v <- calls[[task_names[t]]][id]
          if (is.na(v)) stop(sprintf("glow_scan.wdl, task %s: placeholder ~{%s} is neither a task input of the call nor a known task variable", task_names[t], id))
          w <- regmatches(v, gregexpr("[A-Za-z_][A-Za-z0-9_]*", v))[[1L]]
          w <- c(w[w %in% inputs$name], unname(DERIVED[w[w %in% names(DERIVED)]]))
          res <- c(res, if (length(w)) w else v)
        }
        literal <- !length(ids)
        maps[[length(maps) + 1L]] <- data.frame(task = task_names[t], field = f,
                                                inputs = paste(unique(res), collapse = ", "),
                                                literal = literal, stringsAsFactors = FALSE)
      }
    }
  }
  maps <- do.call(rbind, maps)
  list(inputs = inputs, maps = maps, calls = calls)
}

# The workflow column of a base setting: the inputs that set it, "task layout" for a literal,
# "base default" when no task writes it.
wdl_column <- function(field, maps) {
  m <- maps[maps$field == field, ]
  if (!nrow(m)) return("base default")
  ins <- unique(unlist(strsplit(m$inputs[nzchar(m$inputs)], ",\\s*")))
  parts <- character()
  if (length(ins)) parts <- sprintf("`%s`", ins)
  if (any(m$literal)) parts <- c(parts, "task layout")
  paste(parts, collapse = "; ")
}

#################### ASSEMBLY ####################

md_table <- function(df) {
  esc <- function(x) { x <- as.character(x); x[is.na(x)] <- ""; gsub("\n", " ", x) }
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  body <- apply(df, 1L, function(r) paste0("| ", paste(esc(r), collapse = " | "), " |"))
  c(hdr, sep, body)
}
code <- function(x) ifelse(nzchar(x), paste0("`", gsub("`", "'", x), "`"), "")
NOTES <- new.env()                          # long defaults, shown in full as footnotes
short_default <- function(d, name) {
  if (nchar(d) <= 60L) return(code(d))
  assign(name, d, envir = NOTES)
  paste0(code(paste0(substr(d, 1L, 57L), "...")), " [^", name, "]")
}

build <- function() {
  wdl <- parse_wdl(WDL_FILE)
  cohort_names <- parse_assignments(COHORT_FILE)$name
  all_scripts <- template_scripts(c("00-data-prep", "01-training", "02-single-variant", "03-snv-set", "04-summary"))
  all_codes <- lapply(all_scripts, code_of)
  base_names <- character(); shown_names <- character(); all_docs <- list(); out <- character()
  out <- c(out,
    "# GLOWanalyses settings reference",
    "",
    "Generated by `base-config/make-settings-reference.R` from the four base configuration files,",
    "the shipped example (`runs/example/_cohort.R` and its run configurations), the template scripts",
    "and `cloud/terra/glow_scan.wdl`. Do not edit it by hand: regenerate it from the GLOWanalyses root",
    "with `Rscript base-config/make-settings-reference.R`; the image build runs `--check` and fails",
    "when this file is stale.",
    "",
    "A run's `config.R` is `source(base); source(_cohort); <overrides>`, so a field set later wins",
    "(README, \"Run organization\"). The columns:",
    "",
    "- **Default**: the value the base file assigns; \"required\" for a field the base file documents",
    "  but a run must set.",
    "- **Set by**: where the shipped example sets the field: `cohort` (`_cohort.R`, the data), `run`",
    "  (the run's `config.R`, what the run varies), or `base` (the default is kept).",
    "- **Read by**: the template scripts whose code names the field (a word match; the stage",
    "  directory is the table's heading).",
    "- **Workflow**: on Terra, the input of `glow_scan.wdl` that sets the field; \"task layout\" when",
    "  the task writes a fixed value (a path under `/tmp/glow`); \"base default\" when no task writes",
    "  it, so the workflow runs the base value.",
    "- **Description**: the base file's comment.",
    "")
  counts <- list()
  for (st in names(STAGES)) {
    s <- STAGES[[st]]
    base <- parse_assignments(s$file)
    base_names <- c(base_names, base$name)
    docs <- doc_list(s$file)
    run_names <- unique(unlist(lapply(s$runs, function(f) parse_assignments(f)$name)))
    scripts <- template_scripts(s$dirs); codes <- lapply(scripts, code_of)
    # required fields: documented, not assigned in the base, set by the example or read by the stage
    req <- setdiff(names(docs), base$name)
    req <- req[req %in% c(cohort_names, run_names)]            # a required field is one the example's layers set
    rows <- rbind(
      data.frame(name = base$name, default = base$default, description = base$description, required = FALSE, stringsAsFactors = FALSE),
      data.frame(name = req, default = "", description = unname(docs[req]), required = TRUE, stringsAsFactors = FALSE))
    rows$description <- ifelse(nzchar(rows$description), rows$description,
                               ifelse(rows$name %in% names(docs), unname(docs[rows$name]), ""))
    set_by <- vapply(rows$name, function(n) {
      v <- c(if (n %in% cohort_names) "cohort", if (n %in% run_names) "run"); if (!length(v)) "base" else paste(v, collapse = ", ") }, character(1))
    tbl <- data.frame(
      Setting     = code(rows$name),
      Default     = ifelse(rows$required, "required", vapply(seq_len(nrow(rows)), function(i) short_default(rows$default[i], rows$name[i]), character(1))),
      `Set by`    = set_by,
      `Read by`   = vapply(rows$name, function(n) fmt_readers(readers_of(n, scripts, codes)), character(1)),
      Workflow    = if (st %in% c("data_prep", "training")) "not run by the workflow"
                    else vapply(rows$name, function(n) wdl_column(n, wdl$maps), character(1)),
      Description = rows$description,
      check.names = FALSE, stringsAsFactors = FALSE)
    counts[[st]] <- c(base = nrow(base), required = length(req), nodesc = sum(!nzchar(rows$description)))
    shown_names <- c(shown_names, rows$name); for (n in names(docs)) if (is.null(all_docs[[n]])) all_docs[[n]] <- docs[[n]]
    out <- c(out, sprintf("## %s (`%s`)", s$title, s$file), "",
             sprintf("Templates: %s.", paste(sprintf("`%s/`", s$dirs), collapse = " and ")), "",
             md_table(tbl), "")
  }
  # in-code optional fields
  inc <- incode_fields(all_scripts, all_codes)
  inc <- inc[!inc$name %in% shown_names, ]
  inc <- inc[order(inc$name, inc$script), ]
  inc <- do.call(rbind, lapply(split(inc, inc$name), function(d) data.frame(
    name = d$name[1L], script = paste(d$script, collapse = ", "),
    default = paste(unique(d$default), collapse = "; "),
    description = if (is.null(all_docs[[d$name[1L]]])) "" else all_docs[[d$name[1L]]], stringsAsFactors = FALSE)))
  out <- c(out, "## Settings read with an in-code default", "",
           "These fields have no line in a base file. A template reads them with `g0()` or `get0()` and",
           "uses the default shown when the configuration does not set them; a run sets them in its",
           "`config.R` like any other field.", "",
           md_table(data.frame(Setting = code(inc$name), Scripts = inc$script, Default = code(inc$default),
                               Description = inc$description, check.names = FALSE, stringsAsFactors = FALSE)), "")
  # workflow inputs
  wi <- wdl$inputs
  sets <- vapply(wi$name, function(n) {
    f <- unique(wdl$maps$field[vapply(strsplit(wdl$maps$inputs, ",\\s*"), function(v) n %in% v, logical(1))])
    if (length(f)) paste(code(f), collapse = ", ")
    else if (grepl("^(cpu|mem_gb_|disk_|preemptible|docker_image)", n)) "task runtime"
    else if (grepl("^(Array\\[)?File", wi$type[wi$name == n])) "file input, linked into the task layout"
    else "workflow control" }, character(1))
  out <- c(out, "## Workflow inputs (`cloud/terra/glow_scan.wdl`)", "",
           "Every input of the workflow, in the order of its `input` block, with the type and default the",
           "WDL declares, the configuration field it sets (through the task's `config.R`), and the",
           "description Terra shows in its input table. A field that no input sets runs with its base",
           "default on Terra (the \"Workflow\" column above).", "",
           md_table(data.frame(Input = code(wi$name), Type = code(wi$type), Default = code(wi$default),
                               Sets = sets, Description = wi$description, check.names = FALSE, stringsAsFactors = FALSE)), "")
  notes <- sort(ls(NOTES))
  if (length(notes)) {
    out <- c(out, "## Long defaults", "")
    for (n in notes) out <- c(out, sprintf("[^%s]: `%s`", n, gsub("`", "'", get(n, envir = NOTES))), "")
  }
  list(text = unname(out), counts = counts, n_inputs = nrow(wi), n_incode = nrow(inc), maps = wdl$maps, base_names = base_names)
}

#################### MODES ####################

b <- build()
if (mode == "write") {
  writeLines(b$text, OUT_FILE)
  cat(sprintf("wrote %s: %d base settings, %d required, %d in-code, %d workflow inputs\n", OUT_FILE,
              sum(vapply(b$counts, `[[`, numeric(1), "base")), sum(vapply(b$counts, `[[`, numeric(1), "required")), b$n_incode, b$n_inputs))
} else if (mode == "check") {
  if (!file.exists(OUT_FILE)) { cat("STALE: ", OUT_FILE, " does not exist; run the generator\n"); quit(status = 1L) }
  cur <- readLines(OUT_FILE, warn = FALSE)
  if (identical(unname(cur), b$text)) { cat("settings-reference.md is current\n") } else {
    n <- max(length(cur), length(b$text)); d <- which(c(cur, rep("", n - length(cur))) != c(b$text, rep("", n - length(b$text))))
    cat(sprintf("STALE: %s differs from the generated text at %d line(s); first: line %d\n  committed: %s\n  generated: %s\n",
                OUT_FILE, length(d), d[1L], if (d[1L] <= length(cur)) cur[d[1L]] else "<missing>", if (d[1L] <= length(b$text)) b$text[d[1L]] else "<missing>"))
    quit(status = 1L)
  }
} else {
  # invariants, against counts computed without the parser
  fail <- character()
  for (st in names(STAGES)) {
    f <- STAGES[[st]]$file
    n_regex <- sum(grepl("^[A-Za-z_][A-Za-z0-9_.]*\\s*<-", readLines(f, warn = FALSE)))
    if (n_regex != b$counts[[st]][["base"]]) fail <- c(fail, sprintf("%s: %d assignments by regex, %d parsed", f, n_regex, b$counts[[st]][["base"]]))
  }
  src <- readLines(WDL_FILE, warn = FALSE); wf <- grep("^workflow GLOWScan \\{", src)
  inp <- src[block_lines(src, "^  input \\{", "^  \\}$", from = wf)]
  n_inputs_regex <- sum(grepl("^\\s+[A-Za-z]", inp))
  if (n_inputs_regex != b$n_inputs) fail <- c(fail, sprintf("workflow inputs: %d by regex, %d parsed", n_inputs_regex, b$n_inputs))
  heads <- c(grep("^## ", b$text), length(b$text) + 1L)
  for (st in names(STAGES)) {                                 # each base setting once within its stage's table
    h <- grep(paste0("^## .*\\(`", STAGES[[st]]$file, "`\\)"), b$text)
    sect <- b$text[h:(heads[heads > h][1L] - 1L)]
    for (n in parse_assignments(STAGES[[st]]$file)$name)
      if (sum(grepl(paste0("^\\| `", n, "` \\|"), sect)) != 1L) fail <- c(fail, sprintf("%s: setting %s does not appear exactly once in its table", STAGES[[st]]$file, n))
  }
  wi_txt <- b$text[grep("^## Workflow inputs", b$text):length(b$text)]
  if (sum(grepl("^\\| `[A-Za-z_]+` \\| `", wi_txt)) != b$n_inputs) fail <- c(fail, "the workflow-inputs table does not hold one row per input")
  nodesc <- sum(vapply(b$counts, `[[`, numeric(1), "nodesc"))
  cat(sprintf("self-test: %d base settings, %d in-code settings, %d workflow inputs, %d settings without a description\n",
              sum(vapply(b$counts, `[[`, numeric(1), "base")), b$n_incode, b$n_inputs, nodesc))
  if (length(fail)) { cat(paste0("FAIL: ", fail, "\n"), sep = ""); quit(status = 1L) }
  cat("self-test: PASS\n")
}
