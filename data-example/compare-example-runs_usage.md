
# `compare-example-runs.R`: compare two runs of the shipped example

## What it does

The script compares two runs of the example, file by file and cell by cell. It answers one
question: did a change of the workflow leave the results of the example as they were? The change
can be a new release, the YAML form of the configuration, or a new way of running the stages.

Each run is a GLOWanalyses directory in which the example was run. A run has one of two forms.
`R` is the form of release 0.3.6, with one `config.R` per stage and the outputs beside it. `yaml`
is the form of one YAML configuration, with the outputs on the derived layout under
`runs/example/`. The script maps the folders of the two forms onto each other, so the two runs
may have different forms.

What it compares, and how:

- **A table** (`.csv`, `.tsv`): the same columns, the same rows, and every cell identical as
  text. A missing cell has to be missing in both runs.
- **An R object that is a result** (`.rds`): identical. An object that holds a function is never
  identical to a second copy of itself, so it passes when `all.equal()` with tolerance 0 finds
  no difference, and the report says so.
- **A GDS file**: the same nodes, each with identical data.
- **The two run records of an SNV-set run** (`01-config_used.rds`, `01-input_checksums.rds`):
  field by field, and input by input.

Some values describe the run and not the data: a creation time, the name of the run, a package
version, a folder on the machine. They cannot be equal between two runs. The script lists them
at its top (`RUN_KEYS`, `RUN_COLUMNS`, `OBJECT_RUN_FIELDS`, `GDS_RUN_ATTRIBUTES`), leaves them out
of the comparison, prints both values, and checks the keys on run B. A run name has to be the
name of run B. A path has to exist. The GDS folder and the phenotype bundle have to be those of
run B, so a record that points at another run's inputs fails. Logs, provenance files,
the completion records and the resolved configuration of a YAML run, and plots are counted and
not compared. A plot is a drawing of a table that is compared. The LD scores of a single
chromosome exist only in a run that computed them one chromosome at a time. They are compared
when both runs have them. A run that lacks them passes only when both runs hold the combined LD
scores and these were compared as equal.

## What a PASS needs

No difference is not enough for a PASS, because two runs that both stopped early do not differ.

- **The required results.** Each group has to hold the results that a finished run of the
  example produces, in both runs (`REQUIRED` at the top of the script). A required result that
  was not compared is a failure with its own row in the report. Two kinds of result are
  required by what the run states itself: every chunk of the table that the preparation of a
  window scan wrote, and every model that the PI training says it trained. So two runs that
  both lack the second chunk, or four of five models, do not pass.
- **Checksum records that are true.** A checksum record first has to list the inputs that the
  settings record beside it names: the phenotype bundle, the B model, the files of the PI model
  folder and the GDS file of each chromosome, and no other file. It then has to state the
  checksum of each file it names, in its own run. Only then are the two records compared.
- **The GLOWr version of a file.** GLOWr writes its version into the B model and the annotated
  GDS files. Two runs with two GLOWr versions differ there. The value is left out only when it
  is the version that the run itself recorded (`01-run_info.tsv`), so a file from another run
  does not pass as this run's.

## How to run it

```bash
Rscript data-example/compare-example-runs.R \
  --a <GLOWanalyses dir of run A> --a-form R \
  --b <GLOWanalyses dir of run B> --b-form yaml \
  --out comparison.md
```

`--groups` chooses what is compared. The default is `single_variant,gene,window,coding`, the
four runs of the example. Add `training,data_prep` when both runs rebuilt the steps before them:

```bash
Rscript data-example/compare-example-runs.R --a <dir A> --a-form R --b <dir B> --b-form yaml \
  --groups single_variant,gene,window,coding,training,data_prep --out comparison.md
```

Two runs on one machine give the same text in every cell, and the default asks for that. Two
runs on different machines can differ in the last digits of a number. `--tolerance 1e-8` then
compares numbers by their relative difference and everything else as text, and the report says
so. The rule is the same in a table and in an R object: each number is compared with its
partner, by |a - b| / max(|a|, |b|), and never as a mean over a vector. `--a-root` and `--b-root` give the output root of a YAML run when it is not `runs/` in its
GLOWanalyses directory.

The script needs the R packages `data.table`, `gdsfmt` (for GDS files) and `GLOWr` (to compare
two variant filters).

## What it prints

The report is Markdown: a summary, one row per file with its verdict, and the values that
describe the run. The last line on the screen is `Verdict for the comparison: PASS` or
`Verdict for the comparison: FAIL`, and the exit status is 1 on a failure. A file that exists in
one run only is a failure, and so is a run without any result file.

## Inputs that two runs built separately

An SNV-set run records a checksum of each input: the phenotype bundle, the B model, the PI
models and the GDS files. When two runs built these inputs themselves, a checksum can differ
although the content is the same, because an input carries a label of its run, a time or a
folder. The script then looks for that input in the same report: a row whose two files have
exactly the checksums that the two records state. The checksum record passes only when every
such input was compared and differs in run values only. That needs the
groups `training` and `data_prep`, so a comparison of two rebuilt runs without them fails and
says which groups to add.
