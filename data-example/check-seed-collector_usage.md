# check-seed-collector.sh

`data-example/check-seed-collector.sh` is the regression check of
`data-example/write-seed-settings.R --collect`, the tool that copies the rebuilt pre-built
files of the example into `data-example/cohort/`. It is for a maintainer who changes that
tool, the completion records or the seeding of the example.

## What it checks

The collector may publish a file only when it is the result of a finished unit of
`example.yml` as it is, and the published files and `data-example/seed-settings.yml` have to
be one set. The script disturbs a collection in the ways a second process could, one at a
time, and expects that nothing is collected and that the tracked files keep their bytes:

| Disturbance | What has to happen |
| --- | --- |
| `example.yml` is edited while the files are copied | the collector stops |
| a collected unit is run again meanwhile (its record is gone) | the collector stops, and the unfinished bytes are not published |
| a raw input of the example is changed meanwhile | the collector stops |
| a part before the collected files is run again meanwhile | the collector stops |
| `--collect` together with `--check`, and an unknown argument | the script stops and changes nothing |
| the collection is interrupted after its first file | the marker stays, `--check` fails, and a run of the example does not use the files |
| a pre-built file is replaced by another file | `--check` fails, and a run of the example names the file |

An undisturbed collection then has to succeed, remove the marker and give the same bytes.

## Usage

Run from the GLOWanalyses directory, with the R packages of the README on the library path.

```bash
bash data-example/check-seed-collector.sh /scratch/collector-check   # a folder that is empty or does not exist
```

The script copies the tree into `<work folder>/ga` and works only there. It first runs
`data-example/build-intermediates.sh` in the copy (the data preparation and the training of
the example, a few minutes). The disturbances are made by `data-example/seed-collector-probe.R`.

## Output

One line per check, then `Checks: <n>, failed: <m>.` and
`Verdict for the seed collector: PASS` or `FAIL` (exit status 1). The log of every step is in
the work folder.
