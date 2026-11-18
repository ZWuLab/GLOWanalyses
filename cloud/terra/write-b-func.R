#!/usr/bin/env Rscript
# ==============================================================================
# File: cloud/terra/write-b-func.R
# Purpose: Put the text of a B function into the YAML configuration that a task of
#   cloud/terra/glow_scan.wdl has written. The workflow passes the function as a
#   string input. Written into the configuration by the workflow engine, between the
#   lines of a YAML block, the text depended on how the engine indents what it
#   inserts: a function over several lines got four more blanks at the start of
#   every line after the first, which changes a text constant that runs over a line
#   break, and with it the function. So the engine writes the text to a file of its
#   own (write_lines), and this script writes the setting with the YAML library that
#   the loader reads with. The text that the loader gets is the text of the input.
#
# Usage (from the GLOWanalyses directory, inside a task):
#   Rscript cloud/terra/write-b-func.R <config.yml> <file with the text of the function>
#
#   <config.yml>  the task's configuration. It has to hold the line `  b_func: ~`
#                 once, in its section `snv_set`. The line is replaced.
#   <file>        the text, as write_lines([b_func]) of the workflow wrote it. One
#                 line break at its end belongs to the file and not to the text. A
#                 line end of two characters (CR LF, from a text typed on Windows)
#                 is written as one line break, which is how R reads a script; R
#                 does not parse a function text with CR LF line ends otherwise.
#
# The script stops when the line is missing or there twice, when the text is empty,
# and when the configuration read back does not hold the text byte for byte.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript cloud/terra/write-b-func.R <config.yml> <file with the text of the function>", call. = FALSE)
config <- args[1L]; text_file <- args[2L]
for (f in c(config, text_file)) if (!file.exists(f)) stop("write-b-func.R: file not found: ", f, call. = FALSE)

# the text, as bytes: no line is split, joined or re-encoded on the way
text <- readChar(text_file, file.size(text_file), useBytes = TRUE)
if (length(text) != 1L || is.na(text)) text <- ""
text <- gsub("\r\n", "\n", text, fixed = TRUE)     # a line end of CR LF is a line break, as when R reads a script
text <- sub("\n$", "", text)                       # the line break that ends the file
if (!nzchar(trimws(text))) stop("write-b-func.R: the text of the B function is empty.", call. = FALSE)
if (!validUTF8(text)) stop("write-b-func.R: the text of the B function is not valid UTF-8.", call. = FALSE)
Encoding(text) <- "UTF-8"

lines <- readLines(config, warn = FALSE)
at <- which(lines == "  b_func: ~")
if (length(at) != 1L)
  stop("write-b-func.R: ", config, " has to hold the line `  b_func: ~` once (found ", length(at), " time(s)).", call. = FALSE)
section <- utils::tail(grep("^[A-Za-z_]+:\\s*$", lines[seq_len(at)], value = TRUE), 1L)
if (!identical(section, "snv_set:"))
  stop("write-b-func.R: the line `  b_func: ~` of ", config, " is not in the section `snv_set`.", call. = FALSE)

# The setting as the YAML library writes it, moved under its section by two blanks on every
# line. A text over several lines becomes a block or a quoted text with its line breaks marked,
# whichever the library chooses: it is the library that reads it back.
encoded <- strsplit(yaml::as.yaml(list(b_func = text)), "\n", fixed = TRUE)[[1L]]
out <- c(lines[seq_len(at - 1L)], paste0("  ", encoded), lines[seq_len(length(lines) - at) + at])
tmp <- paste0(config, ".tmp-", Sys.getpid())
writeLines(out, tmp, useBytes = TRUE)

# read back: the same text, and every other line as it was
back <- yaml::read_yaml(tmp)$snv_set$b_func
if (!identical(enc2utf8(back), enc2utf8(text))) {
  unlink(tmp, expand = FALSE)
  stop("write-b-func.R: the configuration read back does not hold the text of the B function as it was given.", call. = FALSE)
}
if (!file.rename(tmp, config)) { unlink(tmp, expand = FALSE); stop("write-b-func.R: could not write ", config, call. = FALSE) }
cat(sprintf("write-b-func.R: the B function (%d line(s), %d character(s)) is written into %s\n",
            length(strsplit(text, "\n", fixed = TRUE)[[1L]]), nchar(text), config))
