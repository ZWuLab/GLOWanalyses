# GLOW in the cloud: the Docker image and the Terra workflow

GLOW runs on Terra, and on any machine that runs WDL workflows in containers, from the same
scripts as the SLURM route. This folder holds the two pieces that make that possible. The first
is the recipe of the official Docker image, which carries R, the GLOW packages and this
GLOWanalyses tree at matching versions. The second is the workflow `terra/glow_scan.wdl`, which
runs the scan stages of GLOWanalyses as thin tasks over those scripts: the single-variant scan
with its calibration, and the gene, window, coding or non-coding SNV-set scan with its
aggregation. Both routes call the shipped analysis templates, so they write the same tables.

What has been verified so far. On the bundled synthetic two-chromosome cohort, run with
miniwdl under Docker and under Apptainer, the workflow's single-variant, gene, window and
coding tables equal `run-example.sh`'s, and a gene run calibrated by the single-variant run
carries its inflation factor, whether it takes the two files or the one-file report. The checks
are `terra/check-wdl-example.R`, run by `terra/run-example-wdl.sh`. A run on Terra itself and
the non-coding masks are not yet part of that record.

What you need: per-chromosome annotated GDS files (the aGDS that `00-data-prep` builds), the
phenotype bundle, a PI ensemble as a `.tar.gz` of `model_*.rds` files, and a B model (or a
closed-form B function). The data-preparation and training stages are not in the workflow yet.
Run them on a cluster, or ask the authors for trained models.

## Running on Terra

1. **Import the workflow.** In a Terra workspace, add `cloud/terra/glow_scan.wdl` as a workflow
   (from Dockstore once the workflow is registered there, or by uploading the file). Every task
   names the image `ghcr.io/zwulab/glowanalyses:<version>`, the one published with this version
   of GLOWanalyses, so nothing is installed at run time.
2. **Put the inputs in a data table, one row per cohort**, with the columns of the table below.
   Then upload one pair of configuration files from `cloud/terra/configs/` on the workflow's
   configuration page. The inputs file binds the workflow inputs to those columns, and the
   outputs file writes each named output back to a column of the same row. There is one pair
   per run type, and each type writes its own columns, so the runs of a cohort do not overwrite
   each other.
3. **Run the single-variant scan first**, with the `GLOWScan_single_variant` pair. It writes the
   raw and the calibrated single-variant tables, the LD scores, the plots and the calibration
   report, and the report lands in the `calibration_report` column.
4. **Run the SNV-set scans**, with the `GLOWScan_gene`, `GLOWScan_window`, `GLOWScan_coding` or
   `GLOWScan_noncoding` pair. Each reads the PI tarball, the B model and the
   `calibration_report` column, so the steps chain without copying paths. The report calibrates
   GLOW's statistics as `03-snv-set/prepare.R` does, and a run without any calibration input is
   uncalibrated, as on the cluster.

**The first run.** Start with one chromosome and `region_type = "gene"` without calibration
inputs. It exercises localization, the fixed paths, the model loading, the shared Stage-1
artifacts, the scan, the aggregation and the plots in a few minutes. A one-chromosome
single-variant run needs `calibration_method = "lambda_gc"`: the default `ldsc_intercept`
regresses across chromosome blocks and stops with one.

**Before a real run, check the defaults.** The variant filter, the window geometry and the
minimum set size follow GLOW's full-spectrum design, not STAARpipeline's rare-variant
conventions. "Defaults to check before a real run" below lists each default and the exact
lines that switch either route to the STAARpipeline settings.

### The data table and the configuration files

The configuration files are in the format Terra's own "Download json" writes: a value that
reads a column of the data table is written `${this.<column>}`, and a literal is a JSON value.
Terra's input configuration documentation describes the `this.` expressions. The inputs files
read these columns:

| Column | Type | What it holds |
|---|---|---|
| `gds_files` | array of files | the per-chromosome annotated GDS files, in the order of `chromosomes` |
| `chromosomes` | array of strings | the chromosome numbers, 1 to 22, parallel to `gds_files` |
| `pheno_bundle` | file | the phenotype bundle `.rds` |
| `pi_models_tar` | file | the PI ensemble as a `.tar.gz` of `model_*.rds` (set-based runs) |
| `b_model_rds` | file | the B model `.rds` (set-based runs) |
| `calibration_report` | file | written by the single-variant run, read by the set-based runs |

The outputs files write to columns named `glow_<run type>_<output>`, for example
`glow_gene_results`, `glow_noncoding_completeness` and `glow_window_loci`, with the one
exception above: the single-variant run's `calibration_report_out` goes to `calibration_report`.
The output names are listed under "Outputs" below. An analysis setting that is not in an inputs
file keeps the workflow's default. To change one, add a line such as
`"GLOWScan.staar_enabled": true` to the file, or set it on the configuration page.

### Inputs

The analysis settings carry the names of the configuration fields in the GLOWanalyses README
("Configuration reference"), so that documentation serves both routes. The file inputs are
localized by the engine and translated into the template paths:

| Workflow input | What the template reads |
|---|---|
| `gds_files`, `chromosomes` | `gds_dir` and `gds_pattern`: the i-th file is linked as `chr<N>.gds` for the i-th chromosome number, so the files may have any basenames; the numbers are 1 to 22, distinct, and the two arrays have the same length |
| `pheno_bundle` | `pheno_path`: the `.rds` list with `sample_id`, `Y` (the trait), `X` (the covariate matrix) and `trait` (`binary` or `continuous`), as `00-data-prep/assemble-pheno-covar.R` writes it |
| `pi_models_tar` | `pi_model_dir`: the `model_*.rds` files of the PI ensemble, at the archive root or inside one directory; required for every set-based run |
| `b_model_rds`, `b_func` | `b_model_path`, or `b_func` when no model is given |
| `calibration_report` | `calibration$report`: the one-file calibration, the `single_variant_calibration.tsv` of a single-variant run on the same phenotype bundle. Its factor and method are taken as estimated there, and Stage 1 refuses a report whose `pheno_bundle_md5` is not this run's bundle, because the factor is specific to the samples, outcome and covariates |
| `single_variant_scan`, `ld_scores`, `calibration_method`, `inflation_factor` | the two-file alternative (`marginal_all.csv` and `ld_scores.rds` of that run, the factor re-estimated by `calibration_method`), or a manual `inflation_factor`. One of the three forms, or none: without any of them the run is uncalibrated |
| `cpu`, `n_cores` | `cpu` is the cores requested; the workers are `n_cores`, which every scan task takes from the vCPUs of the machine it runs on unless the input is set (see "Sizing a run") |
| `use_spa` | `use_spa`: `auto` (SPA for a binary trait), `spa` or `standard` |
| `staar_enabled` | the embedded STAAR comparator on GLOW's own variant sets; coding runs with it need a binary trait. The native STAARpipeline stage of the SLURM route is not part of the workflow |

`terra/glow_scan.example.*.inputs.json` are complete input sets for the bundled synthetic cohort,
for miniwdl. Their paths are relative to a GLOWanalyses checkout. Every input has a description
in the workflow's `parameter_meta` block, which Terra shows in its input table.

### Defaults to check before a real run

The defaults follow GLOW's design, which allows for common variants, rather than limiting the
analysis to rare variants as a rare-variant pipeline does (for example, the typical STAARpipeline
setting). A run that keeps them answers a different question than a STAARpipeline run of the same
data, so check them before a real run. You can certainly reset them to match your analysis needs:

- `rare_maf_cutoff = 0.5` keeps the whole frequency spectrum by default. For a sequencing-style
  rare-variant analysis, set `0.01`, STAARpipeline's cutoff.
- `window_size = 100000` and `step_size = 10000`: 100 kb windows every 10 kb. STAARpipeline's
  sliding windows are 2 kb every 1 kb. Set `2000` and `1000` to follow them.
- `staar_enabled = false`. The STAAR comparator adds STAAR's test to every cell, so a run with it
  takes longer. Turn it on when the STAAR columns are wanted beside GLOW's.
- `min_variants = 1`: a one-variant set is tested as its single variant and flagged
  `is_single`. STAARpipeline tests sets of two or more variants (its `rv_num_cutoff`). Set `2`
  to follow it.

The defaults are the same on both routes, because the workflow mirrors the SLURM route's base
configuration, `base-config/snv_set_base.R`. Stage 1 freezes the settings into its snapshot:
every Terra workflow run applies them in its first task, and on the cluster a change made after
Stage 1 needs `03-snv-set/prepare.R` again. A finished run records the values it used in
`01-run_info.tsv` (the `run_info` output).

**On Terra.** Add the settings you change to the inputs file of the run type before you upload
it on the workflow's configuration page (each line inside the braces, with a comma after every
line but the last), or type the value into the input's field on that page. The STAARpipeline
settings, all five:

```json
"GLOWScan.rare_maf_cutoff": 0.01,
"GLOWScan.min_variants": 2,
"GLOWScan.window_size": 2000,
"GLOWScan.step_size": 1000,
"GLOWScan.staar_enabled": true
```

`window_size` and `step_size` are read by window runs only (`GLOWScan_window.inputs.json`). The
other three apply to every set-based pair.

**On the cluster.** Put the same settings in the run's own configuration,
`runs/<name>/snv-set/config.R` (`snv-set-window/config.R` for a window run), after its two
`source()` lines. Leave `base-config/snv_set_base.R` as shipped (the GLOWanalyses README, "Run
organization"). `filter_spec` is one list, so restate all four of its fields:

```r
filter_spec   <- list(rare_maf_cutoff = 0.01, variant_type = "SNV",
                      min_mac = 1L, min_variants = 2L)
staar_enabled <- TRUE      # the STAAR package must be installed
window_size   <- 2000L     # window runs only
step_size     <- 1000L
```

`runs/example/snv-set-window/config.R` is a worked example of such an override block.

### Outputs

Each output is one table or one array of per-chromosome files, under the template's own file
name. The configuration files map every output of a run type to a column. A run of one type
leaves the other types' outputs null.

Every run:

| Output | Content |
|---|---|
| `per_chrom_logs` | the template's log of each per-chromosome task |
| `per_chrom_resources` | the resource monitor's lines of each per-chromosome task (see "Sizing a run") |
| `run_tag` | `region_type`, or `region_type_<run_label>` when a label is given; the file names do not change with it |

A single-variant run:

| Output | File |
|---|---|
| `single_variant_results` | `marginal_all.csv`, the genome-wide raw table |
| `single_variant_calibrated` | `single_variant_all_calibrated.csv`, the calibrated p-values |
| `calibration_report_out` | `single_variant_calibration.tsv`: the factor at full precision, its method, the regression's size, the raw and calibrated lambdas, the phenotype bundle's checksum and the versions; the one input a set-based run's calibration takes |
| `single_variant_calibration` | `single_variant_calibration.txt`, the same facts as text |
| `ld_scores_out` | `ld_scores.rds` |
| `single_variant_plots` | `marginal_plots.pdf`, the Manhattan and QQ plots |
| `calibrate_log`, `calibrate_resources` | the combining task's log and monitor lines |

A set-based run:

| Output | File |
|---|---|
| `run_info` | `01-run_info.tsv`: the run as one key/value table (type, chromosomes, cohort size and trait, filter, inflation factor with its source, settings, models, versions) |
| `shared_artifacts`, `prepare_log`, `prepare_resources` | the Stage-1 snapshot and models as a tarball, with the task's log and monitor lines |
| `per_chrom_results` | the per-chromosome tables (`glow_chr<N>.csv`, `scan_chunk_<k>.csv`, `glow_coding_chr<N>.csv` or `glow_noncoding_chr<N>.csv`) |
| `per_chrom_errors` | the `*_errors.csv` sidecars: the cells that errored, with the message per cell (a non-coding unit writes one per chromosome, with zero rows when clean; the other types only when a cell errored) |
| `per_chrom_summaries` | non-coding only: the `*_summary.csv` per-mask tallies |
| `combined_results` | the genome-wide table: `glow_results_all.csv` (gene), `scan_results_all.csv` (window), `glow_coding_all.csv` (coding) or `glow_noncoding_all.csv` (noncoding) |
| `top_hits` | gene: `glow_top_hits.csv` |
| `sig_windows`, `loci` | window: `scan_sig_windows.csv`, `scan_loci.csv` |
| `coding_merged`, `coding_sig_by_category`, `coding_sig_cells` | coding: `coding_merged.csv`, `coding_sig_by_category.csv`, `coding_sig_by_category_cells.csv` |
| `noncoding_merged`, `noncoding_sig_by_category`, `noncoding_sig_cells`, `noncoding_thresholds`, `noncoding_completeness` | noncoding: the merged table, the two significance tables, the Bonferroni cutoffs per mask and overall, and the completeness table |
| `lambda_by_category` | coding and noncoding: `lambda_by_category.csv` |
| `plots` | the PDFs of `04-summary/plots.R` (the Manhattan and the QQ plots) |
| `aggregated_tar` | the whole `04-summary` folder as one archive |
| `aggregate_log`, `aggregate_resources` | the aggregation task's log and monitor lines |

### Sizing a run: workers, memory and the resource monitor

**Workers.** Every scan task uses the vCPUs of the machine it runs on as its forked workers,
not the `cpu` request. A cloud platform can give more vCPUs than requested when the memory
request needs them (on Google Cloud a machine has at most 6.5 GB of memory per vCPU), and those
cores are paid for whether they work or not. The task's log opens with a line
`workers: N (vCPUs: N; cpu requested: cpu; memory: M GB)`. Set the `n_cores` input to use fewer
workers than vCPUs, for example when the heaviest genes need more memory per worker.

**Memory.** Size `mem_gb_scan` per vCPU: the per-chromosome tasks run one region per worker at
a time, so their memory grows with the worker count. The templates' size-aware scheduler runs
the heaviest regions first, so a chromosome's wall time is the larger of two numbers, the sum
of its regions' times over the workers and its heaviest region alone. The first task of a
set-based run localizes every chromosome's GDS file, because Stage 1 reads them all. Each
per-chromosome task localizes one. A task's disk is sized from all its localized inputs
(`disk_multiplier`, `disk_extra_gb`).

**The monitor.** Every task runs the monitor of `terra/task-monitor.sh` in the background. It
samples the memory in use every 15 seconds and writes, every 10 minutes and once more at exit,
one line to the task's `resources.txt` (the `*_resources` outputs) and to its standard output:

```
[monitor 2026-10-27T02:36:22Z] memory in use 30.0 GB, peak 30.0 GB, load average 19.8 19.6 18.9, vCPUs 20
[monitor exit 2026-10-27T02:36:23Z] peak memory in use 30.0 GB (this task's processes), machine 128.0 GB, vCPUs 20, exit status 0
```

The memory in use is the memory of the task's own processes, the proportional set size of the
task shell and all its descendants, so the pages that forked workers share are counted once. It
needs no cgroup, which matters on a cluster, where the container shows the whole machine's
cgroup. The peak is the largest sample, so a task that finishes within seconds reports the
memory before its work began. The exit line is written when the command fails too, so an
out-of-memory stop leaves its last readings in the task log. The peak of a representative
chromosome, with a margin, is the number to request for the whole genome. On a machine of its own, as on Terra, a load average near
the vCPU count says every worker is busy.

**Costs.** The hours are the user's cloud bill. The GLOWanalyses README ("Cores and large
genes") has the only published timings, for chromosome 22 of the 1000 Genomes Project. The
monitor lines of your own first chromosome are the sizing for the rest.

### What Stage 1 checks before any scan

A wrong input stops the run in the first task, in minutes, rather than in every chromosome's
scan. The tasks validate their inputs first: parallel, distinct chromosomes 1 to 22; the PI
tarball of a set-based run; `use_spa`; `calibration_method`; one calibration form at most, and
`ldsc_intercept` with `ld_scores`. Then `03-snv-set/prepare.R` checks the data. It counts, per
chromosome, the variants whose `annotation/filter` value is missing (a VCF `FILTER` of `.`), which
count as not passing, and it stops when a chromosome has no QC value at all. It checks the
phenotype bundle against the calibration report. For a non-coding run it checks the five FAVOR
fields the masks read on every chromosome and the two packages the promoter masks need. A cell
that errors during a scan is listed in `per_chrom_errors` with its message, and the scan goes on.

## Running anywhere else with miniwdl

The same workflow runs on a workstation or a cluster node with
[miniwdl](https://miniwdl.readthedocs.io) and either Docker or Apptainer:

```bash
pip install miniwdl                       # into a virtual environment of its own
cd /path/to/GLOWanalyses
# Docker: nothing else to set. Apptainer (clusters), where Docker is absent:
export MINIWDL__SCHEDULER__CONTAINER_BACKEND=singularity
export MINIWDL__SINGULARITY__EXE='["apptainer"]'
export MINIWDL__SINGULARITY__RUN_OPTIONS='["--containall","--no-mount","hostfs"]'   # no --fakeroot
miniwdl run cloud/terra/glow_scan.wdl -i cloud/terra/glow_scan.example.single_variant.inputs.json -d runs/example-wdl/single-variant/.
```

`bash cloud/terra/run-example-wdl.sh` runs the synthetic example through the workflow seven
times: single-variant, gene, gene calibrated by the single-variant run's two files, gene
calibrated by its report, window, coding, and coding with the STAAR comparator. It then compares
the tables with the ones `run-example.sh` wrote in the same tree and checks the named outputs on
disk. `GLOW_WDL_TOLERANCE=1e-12` relaxes the comparison of the two full-precision artifacts (the
LD scores and the calibrated single-variant table) when the reference was computed in another
environment. The p-value tables stay exact. Under Apptainer only the task's working directory and
`/tmp` are writable. The workflow writes nowhere else.

## The image

Two layers, both built by `.github/workflows/docker-image.yaml` at every release tag and
published to the GitHub Container Registry after the example and the workflow comparison pass
on the freshly built image:

- `ghcr.io/zwulab/glow-base:<date>` (`docker/base/Dockerfile`): `rocker/r-ver:4.3.3` with the
  dated CRAN snapshot that image pins (p3m.dev, April 2024), Bioconductor 3.18, and STAAR,
  MultiSTAAR, SCANG and STAARpipeline at pinned commits. The exact commits are the build
  arguments at the top of that Dockerfile. The package manifest is at
  `/opt/glow-base-manifest.tsv` inside the image. The base is rebuilt when `GLOW_BASE_TAG` in
  `docker/versions.env` changes. A hash of the recipe is a label of the image, and the build
  refuses a published base whose label differs from the tree's recipe.
- `ghcr.io/zwulab/glowanalyses:<version>` (`docker/Dockerfile`): the base plus GFisher, GLOWr
  and GLOWpipeline installed from their release-tag archives at the versions named in
  `docker/versions.env`, and this GLOWanalyses tree at `/opt/GLOWanalyses`, the working
  directory, which is where the tasks find the templates and `terra/task-monitor.sh`.
  `/opt/glow-image-manifest.txt` lists what is inside, and every task prints it.

The versions in `docker/versions.env` and the "Requires" line of `NEWS.md` are the same by
construction (the build fails when they disagree). `OPENBLAS_NUM_THREADS=1` is set in the image,
because the scripts parallelize with forked workers. The build context of the official image is
the released tree, which holds no study configuration. `.dockerignore` leaves the
development-only run folders out of a build from a checkout as well. The rocker tag, the system
package repositories and the Bioconductor release are not content-addressed, so the identity of a
tested environment is the image digest and the two manifests, not the recipe alone.

## Not yet in the workflow

The data-preparation stage (VCF or PLINK to GDS, FAVOR annotation, principal components, the
phenotype bundle) and the training stage follow in a later release. The call and the stack of
a cell that errors are not recorded with its message yet. The UK Biobank Research Analysis
Platform is a DNAnexus platform, on which a Terra WDL is not automatically runnable. We have not
tested a route there.

## Credits

The workflow and the image recipes are derived from the Terra workflows written by Alberto
Brusati (Landers lab, UMass Chan Medical School), 2026, who ran GLOW on Terra first and whose
feedback shaped the engine's parallel and large-gene behavior. The one-workflow shape with the
`region_type` switch, the input names and the resource formula are his, and so are the named
outputs with their Terra configuration files, the workers taken from the machine, the resource
monitor, the one-file calibration report and the Stage-1 input checks.
