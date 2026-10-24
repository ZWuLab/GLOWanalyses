# GLOW in the cloud: the Docker image and the Terra workflow

GLOW runs on Terra, and on any machine that runs WDL workflows in containers, from the same
scripts as the SLURM route. This folder holds the two pieces that make that possible. The first
is the recipe of the official Docker image, which carries R, the GLOW packages and this
GLOWanalyses tree at matching versions. The second is the workflow `terra/glow_scan.wdl`, which
runs the scan stages of GLOWanalyses as thin tasks over those scripts: the single-variant scan
with its calibration, and the gene, window, coding or non-coding SNV-set scan with its
aggregation. Both routes call the shipped analysis templates, so they write the same tables.

What has been verified so far. On the bundled synthetic two-chromosome cohort, run with
miniwdl and Docker, the workflow's single-variant, gene, window and coding tables equal
`run-example.sh`'s byte for byte. A gene run calibrated by the single-variant run carries its
inflation factor. The checks are `terra/check-wdl-example.R`, run by
`terra/run-example-wdl.sh`. Runs on Terra itself, with Apptainer on a cluster, and of the
non-coding masks are not yet part of that record.

What you need: per-chromosome annotated GDS files (the aGDS that `00-data-prep` builds), the
phenotype bundle, a PI ensemble as a `.tar.gz` of `model_*.rds` files, and a B model (or a
closed-form B function). The data-preparation and training stages are not in the workflow yet.
Run them on a cluster, or ask the authors for trained models.

## Running on Terra

1. **Import the workflow.** In a Terra workspace, add `cloud/terra/glow_scan.wdl` as a workflow
   (from Dockstore once the workflow is registered there, or by uploading the file). Every task
   names the image `ghcr.io/zwulab/glowanalyses:<version>`, the one published with this version
   of GLOWanalyses, so nothing is installed at run time.
2. **Run the single-variant scan first.** `region_type = "single_variant"` with the GDS files,
   the chromosome numbers and the phenotype bundle. It writes `marginal_all.csv`, `ld_scores.rds`,
   the calibrated table and `single_variant_calibration.txt`.
3. **Run the SNV-set scans.** Set `region_type` to `gene`, `window`, `coding` or `noncoding`,
   and give the PI tarball and the B model. Pass the two files of step 2 as
   `single_variant_scan` and `ld_scores`. They calibrate GLOW's statistics as
   `03-snv-set/prepare.R` does, and without them the run is uncalibrated, as on the cluster.

**The first run.** Start with one chromosome and `region_type = "gene"` without calibration
inputs. It exercises localization, the fixed paths, the model loading, the shared Stage-1
artifacts, the scan, the aggregation and the plots in a few minutes. A one-chromosome
single-variant run needs `calibration_method = "lambda_gc"`: the default `ldsc_intercept`
regresses across chromosome blocks and stops with one.

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
| `single_variant_scan`, `ld_scores`, `calibration_method`, `inflation_factor` | the `calibration` list, or `inflation_factor`. The two files estimate the factor by `calibration_method`. A manual `inflation_factor` is the alternative, and giving both is an error. Without any of them the run is uncalibrated |
| `cpu` | `n_cores`, the forked workers per chromosome |
| `use_spa` | `use_spa`: `auto` (SPA for a binary trait), `spa` or `standard` |
| `staar_enabled` | the embedded STAAR comparator on GLOW's own variant sets; coding runs with it need a binary trait. The native STAARpipeline stage of the SLURM route is not part of the workflow |

`terra/glow_scan.example.*.inputs.json` are complete input sets for the bundled synthetic cohort.
Their paths are relative to a GLOWanalyses checkout, and on Terra the user's bucket files take
their place. Every input has a description in the workflow's `parameter_meta` block, which
Terra shows in its input table.

### Outputs

A single-variant run returns:

- `marginal_all`, `ld_scores_out`, `single_variant_calibrated` and `single_variant_calibration`
  (the factor and the method), the files a set-based run's calibration takes.
- `single_variant_plots` and the per-chromosome `single_variant_logs`.

A set-based run returns:

- `shared_artifacts`, the Stage-1 snapshot and models as a tarball, with `prepare_log`.
- the per-chromosome `scan_results` and `scan_logs`.
- `aggregated`, the tables and plots of `04-summary`, also as `aggregated_tar`, with
  `aggregate_log`.

`run_tag` is `region_type`, or `region_type_<run_label>` when a label is given. The file names
do not change with the label.

**Costs.** The first task of a set-based run localizes every chromosome's GDS file, because
Stage 1 reads them all. Each per-chromosome task localizes one. A task's disk is sized from all
its localized inputs (`disk_multiplier`, `disk_extra_gb`), and `cpu` is the number of forked
workers each chromosome gets (`n_cores`). Memory is the user's to set per task kind. The LD-score
step of a single-variant run holds matrices that grow with the sample count and the variants per
window, so measure a representative chromosome before a whole-genome run. The hours are the
user's cloud bill. The GLOWanalyses README ("Cores and large genes") has the only published
timings, for chromosome 22 of the 1000 Genomes Project.

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

`bash cloud/terra/run-example-wdl.sh` runs the synthetic example through the workflow six times
(single-variant, gene, gene calibrated by the single-variant run, window, coding, and coding
with the STAAR comparator) and then compares the tables with the ones `run-example.sh` wrote in the
same tree. Under Apptainer only the task's working directory and `/tmp` are writable. The
workflow writes nowhere else.

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
  directory. `/opt/glow-image-manifest.txt` lists what is inside, and every task prints it.

The versions in `docker/versions.env` and the "Requires" line of `NEWS.md` are the same by
construction (the build fails when they disagree). `OPENBLAS_NUM_THREADS=1` is set in the image,
because the scripts parallelize with forked workers. The build context of the official image is
the released tree, which holds no study configuration. `.dockerignore` leaves the
development-only run folders out of a build from a checkout as well. The rocker tag, the system
package repositories and the Bioconductor release are not content-addressed, so the identity of a
tested environment is the image digest and the two manifests, not the recipe alone.

## Not yet in the workflow

The data-preparation stage (VCF or PLINK to GDS, FAVOR annotation, principal components, the
phenotype bundle) and the training stage follow in a later release. The UK Biobank Research
Analysis Platform is a DNAnexus platform, on which a Terra WDL is not automatically runnable.
We have not tested a route there.

## Credits

The workflow and the image recipes are derived from the Terra workflows written by Alberto
Brusati (Landers lab, UMass Chan Medical School), 2026, who ran GLOW on Terra first and whose
feedback shaped the engine's parallel and large-gene behavior. The one-workflow shape with the
`region_type` switch, the input names and the resource formula are his.
