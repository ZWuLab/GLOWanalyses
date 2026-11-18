# The synthetic example cohort (`data-example/`)

This folder holds a small, fully synthetic cohort and a small FAVOR database, so that the
GLOWanalyses templates run on a fresh clone without any other data. The example checks the
mechanics of the workflow: every stage writes its tables with the expected columns and at least
one row. Its results carry no statistical meaning.

The example is the dataset `example` of the configuration `example.yml` in the GLOWanalyses
directory. Its results are written to `runs/example/`.

## What is in the folder

This section lists the files that ship and the files that a run writes.

**Shipped** (tracked, about 0.6 MB):

| Path | What |
|---|---|
| `make-example-cohort.R` | the generator of the raw inputs, with a fixed seed |
| `cohort/example_plink/chr{21,22}.{bed,bim,fam}` | the raw genotypes, as PLINK files |
| `cohort/favor-db/` | the small FAVOR database (`FAVORdatabase_chrsplit.csv` and `chr{21,22}_1.csv`) |
| `cohort/pheno.csv` | the phenotype table: `sample_id`, `outcome` (0 or 1), `sex`, `age`, `subpop` |
| `cohort/known-snps.csv` | 20 known variants, for the B model and as the cases of the PI models |
| `cohort/variant_map.csv` | for every variant its gene, its category and its allele frequency |
| `cohort/example_gds_favor/flexible/gds/chr{21,22}_favor.gds` | pre-built: the annotated GDS files |
| `cohort/pheno_covar.rds` | pre-built: the phenotype bundle |
| `cohort/b_model.rds`, `cohort/pi_models/` | pre-built: the B model and five PI models |
| `cohort/pi_models.tar.gz` | the PI models as the one archive that the Terra workflow takes |
| `seed-manifest.yml`, `seed-settings.yml` | which pre-built file stands for which part of the workflow, and the settings the files were built with |
| `build-intermediates.sh`, `write-seed-settings.R` | the maintainer's tools that renew the pre-built files |
| `check-seed-collector.sh`, `seed-collector-probe.R` | the maintainer's check of the second tool: it disturbs a collection in a copy of the tree |
| `check-example.R` | the check that `run-example.sh` runs at its end |
| `compare-example-runs.R`, `check-run-records.sh` | two checks for maintainers, each with a usage note |

**Written by a run** (not tracked): everything under `runs/example/` in the GLOWanalyses
directory.

## The pre-built files

This section says why the example starts at the analysis steps.

The data preparation and the training take the longest and need the most packages. So the
files they produce from the raw inputs are shipped, and the default run of the example starts
at the single-variant step. On its first run, the workflow puts the pre-built files into
`runs/example/` at the places where the two steps would have written them, and gives each one a
record that states the settings it was built with. A run with those settings reads them as
finished results. A run that changes one of those settings, or one of the raw inputs, is told
which step to run first, and that step then produces its own files.

## The cohort

This section describes the synthetic data.

- **Two chromosomes**, 21 and 22, so that the work per chromosome and its combination are
  exercised.
- **300 samples** in two sub-populations with slightly shifted allele frequencies, which gives
  the principal components and the null model something to estimate. About half are cases.
- **374 variants**, 186 on chromosome 21 and 188 on chromosome 22. They lie inside eight real
  genes (`SIK1B`, `CBSL`, `U2AF1L5`, `CRYAA2` on chromosome 21, and `POTEH`, `IL17RA`, `HDHD5`,
  `SLC25A18` on chromosome 22, from the gene table that GLOWr bundles) and between genes. So
  the gene, window and coding scans all form their sets with the default gene table. The
  variants inside genes carry coding categories (plof, missense, disruptive missense,
  synonymous and others). The variants between genes supply the common variants that the
  principal components and the LD-score calibration need.

The annotated GDS files carry the 16 annotation columns that the PI models read and the three
text columns that the coding masks read (`genecode_comprehensive_category`,
`genecode_comprehensive_exonic_category`, `metasvm_pred`), which `example.yml` lists under
`data_prep.favor_features`.

## Running the example

This section gives the commands. They are run from the GLOWanalyses directory.

```bash
bash run-example.sh                  # the single-variant step and the gene, window and coding scans, then the check
REBUILD=1 bash run-example.sh        # also regenerate the cohort and run the data preparation and the training
DRYRUN=1  bash run-example.sh        # show every unit of work and its state, and start none
```

`check-example.R` checks, stage by stage, that the files exist, have the expected columns and
are not empty, that the PI ensemble holds the number of models its training states, and that
every unit of the steps that ran is complete by its record. It exits with a non-zero status on
any failure. It does not judge p-values or AUC values, which mean nothing on synthetic data.

## Renewing the pre-built files

This section is for a maintainer who changes the generator or a default.

```bash
Rscript data-example/make-example-cohort.R      # the raw inputs (the same seed gives the same cohort)
bash    data-example/build-intermediates.sh     # the pre-built files, and seed-settings.yml
```

`build-intermediates.sh` runs the data preparation and the training of `example.yml` with
`--force`, copies their results back into `cohort/`, and writes `seed-settings.yml`. A file is
copied only when a run of `example.yml` as it is now finished the part that wrote it.
`seed-settings.yml` also states the checksum of every pre-built file, and a run of the example
uses the files only when they are that set. While files are put in place, the marker
`seed-collection.incomplete` is there, so a collection that was interrupted is noticed: the
builder is then run again. `write-seed-settings_usage.md` describes that step, and
`bash data-example/check-seed-collector.sh <work folder>` checks it in a copy of the tree
(`check-seed-collector_usage.md`).

