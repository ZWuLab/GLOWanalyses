# write-seed-settings.R

`data-example/write-seed-settings.R` writes `data-example/seed-settings.yml`, the settings
that the example's pre-built files were built with.

## Why the file exists

The example ships files that a full run would produce (the annotated GDS files, the
phenotype bundle, the B model and the PI models), listed in `data-example/seed-manifest.yml`.
On a run of the example, `lib/glow-stage.R` puts them in place and gives each unit a seed
record. A seed record states the settings the files were built with, and a run whose
settings differ is told to produce its own files.

Those settings must not follow a later edit of `example.yml`: the files do not change when
the configuration does. So they are fixed in `seed-settings.yml` at the time the files are
built. For each seeded part the file holds the settings of the part and of every part
before it, and a checksum of the content of their outside files (the raw inputs of the
example). A raw input that is changed in place therefore makes the pre-built files invalid
for a run, and a clone at another path does not.

## Usage

Run from the GLOWanalyses directory.

```bash
Rscript data-example/write-seed-settings.R --collect   # copy the rebuilt files back, then write the file
Rscript data-example/write-seed-settings.R             # write the file from example.yml
Rscript data-example/write-seed-settings.R --check     # compare, write nothing
```

`--collect` is what `data-example/build-intermediates.sh` calls. It takes the folders from
the layout of `example.yml`, and it copies a file of `data-example/seed-manifest.yml` back
into `data-example/cohort/` only when a run of `example.yml` as it is now finished the unit
that wrote it. A unit that is not finished, or that still holds the supplied files, stops
the script before anything is copied. A folder (the PI models) is replaced by the files
that its completion record lists.

Before a file enters `data-example/cohort/`, the script validates every unit a second time,
after it has copied the files into a staging folder, with dependencies that it looks up
anew. A unit, or a part before it, that was run again meanwhile, a raw input that was
changed meanwhile and an `example.yml` that was edited meanwhile all stop the script with
nothing collected. The settings that the file states are those of the inputs that the
validated units were made from.

The file ends with `artifacts`, the checksum of every pre-built file. A run of the example
uses the pre-built files only when they are that set: a file that is missing, a file more
and a file with other bytes stop the run with the name of the file. From the first file
that `--collect` puts in place until the settings file is written, the marker
`data-example/seed-collection.incomplete` is there. A collection that is interrupted leaves
it, and then neither `--check` nor a run of the example accepts the files. Running
`data-example/build-intermediates.sh` to its end removes it.

`--check` prints `Verdict for the seed settings: PASS` when the file states what
`example.yml` and the pre-built files give now and no marker is there, and `FAIL` with exit
status 1 otherwise.

`bash data-example/check-seed-collector.sh <work folder>` is the regression check of the
collector (`check-seed-collector_usage.md`).

## When to run it

Only through `data-example/build-intermediates.sh`, which rebuilds the pre-built files from
`example.yml` and calls this script with `--collect` at its end. Writing the file without
rebuilding the files would claim settings that the files were not built with.

A setting whose default changes in `base-config/settings-registry.yml` changes what
`example.yml` gives, so `--check` then fails until the files are rebuilt. The test suite of
GLOWpipeline runs the same comparison.
