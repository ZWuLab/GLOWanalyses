
# `check-run-records.sh`: check the records of a run on the example

## What it does

A run with a YAML configuration keeps a completion record for every unit of work: what the unit
depended on, which files it produced, and a mark that is new at every run. The records decide
three things. A unit that is finished is skipped. A run that would overwrite results made with
other settings is refused. An analysis that extends another reuses the parts that did not
change.

This script checks those rules on the shipped example. It runs the stage scripts on seven small
configurations that extend `example.yml` and write to a work folder of their own, so
`runs/example/` is not touched. It makes ten groups of checks, 48 checks in all:

- **(a)** A derived analysis that changes the PI features stops at the SNV-set preparation. The
  message names the training step and the setting that differs.
- **(b)** After the PI training has run in the derived analysis, the preparation runs. The B
  model and the phenotype bundle are those of the parent analysis, the PI models are the
  derived analysis's own, and `config-resolved.yml` says so.
- **(c)** A changed `ld_threshold` aimed at the finished folder of the same analysis is refused.
  The message gives the recorded and the current value and the two ways forward: a new
  analysis name, or `--force`. Every result file has the checksum it had before the refusal.
- **(d)** After the forced preparation, the scan units of the finished run are not complete
  any more.
- **(e)** A deleted result file makes its unit run again. The other unit stays complete.
- **(f)** A configuration of another dataset whose raw inputs were left at the example's files
  is stopped, and the message names the setting.
- **(g)** A derived analysis that changes nothing the gene run reads reuses the finished folder
  of its parent, whole: the preparation, a scan unit and the aggregation, and nothing is
  written in the derived analysis.
- **(h)** A parent with a finished preparation and no scan is not a finished folder. The scan
  of an analysis derived from it is not reported as reused from results that do not exist. The
  derived analysis runs its own preparation and scan.
- **(i)** A setting on the command line (`--max-genes`) stops the stage script, names the
  setting to put in the configuration, and leaves no record and no held lock.
- **(j)** A file that is deleted inside a recorded folder (the settings record of the
  preparation) makes the unit run again.

## How to run it

Run it in an R environment that has GLOWr and a GLOWpipeline with the run records
(`glow_unit_begin()`).

```bash
bash data-example/check-run-records.sh              # a new temporary work folder
bash data-example/check-run-records.sh /tmp/records # a work folder of your choice, empty or absent
```

It takes about three minutes. The step (b) trains five small PI models on the example.

## What it prints

Each check prints `[PASS]` or `[FAIL ...]` with its name. A failed check shows the last lines of
its log. The last line is `Verdict for the run records: PASS` or `... FAIL`, and the exit status
is 1 on a failure. The output of every command is kept in the work folder as `NN-<check>.log`.
