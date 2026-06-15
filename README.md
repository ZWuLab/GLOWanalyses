# GLOWanalyses

Runnable, cohort-agnostic **workflow templates** for the GLOW methodology. This
is a *scripts directory, not an R package* (no `DESCRIPTION`/`NAMESPACE`): copy a
run, point it at your data via a config, and run. The statistical methods live in
`GLOWr`; the orchestration + summary live in `GLOWpipeline`; these scripts are
the thin, numbered glue that calls them.

```
00-data-prep/  PLINK→GDS · FAVOR annotation (→aGDS) · PCs · pheno/covar assembly
01-training/   B (effect-size dist.) + PI (pathogenicity ensemble) training
02-single-variant/  per-variant scan · in-sample LD scores · GC calibration
03-snv-set/    Stage 1 (prepare.R) + thin Stage-2 entries (run-gene/window/coding.R)
04-summary/    Stage 3: aggregate.R (region_type-aware) + plots.R (Manhattan/QQ/lambda/drilldown)
base-config/   {data_prep,training,single_variant,snv_set}_base.R - documented defaults (no cohort paths)
data-example/  a self-contained synthetic cohort + generator + self-test (runs out-of-the-box)
runs/<name>/   per-run config.R (sources the base, names the data) + derived outputs/
slurm/         array-job templates (run-gene/window/coding.sh) + submit-throttled.sh
```

All five stages (`00`→`04`) are thin `--config` wrappers over the packaged
`GLOWr` / `GLOWpipeline` methods. The repository ships a **runnable,
self-contained synthetic example** (`data-example/` + `runs/example/`) — run
`bash run-example.sh` from the GLOWanalyses directory and the whole chain
executes on a tiny tracked cohort (see the
[End-to-end example](#end-to-end-example-synthetic-out-of-the-box) walkthrough
below). A `noncoding` SNV-set slot is future work (one more thin entry
+ a `build_scan_regions` dispatcher case, not new machinery).

## Install + run convention

These scripts call the three GLOW packages at runtime; **install them first**
(they are separate repos), then run every command **from the GLOWanalyses
directory** (paths are GLOWanalyses-root-relative):

```r
# in R, once: install the runtime dependencies (GFisher -> GLOWr -> GLOWpipeline)
# e.g. remotes::install_github("ZWuLab/GLOWr") etc., or from local source trees
```

```bash
cd /path/to/GLOWanalyses          # all stage/config paths resolve from here
conda activate r_env              # or any R with the 3 packages installed
```

## The five stages (`00`→`04`)

Each stage is run from the GLOWanalyses directory as `Rscript
<stage>/<script>.R --config <run>/config.R [unit selector]`. The
producer/consumer chain:

| Stage | Scripts | Reads | Writes |
|---|---|---|---|
| **00 data-prep** | `plink-to-gds` · `annotate-favor` · `compute-pcs` · `assemble-pheno-covar` | raw PLINK + a FAVOR DB + a pheno table | the cohort **data-lineage tree** (`<base>_gds`, `_gds_favor` (aGDS), `_pcs`, `_pheno`) |
| **01 training** | `b/{prepare,estimate}` · `pi/{prepare,annotate,train,evaluate}` | a known-SNP table + the FAVOR DB + the aGDS control tree | `b_model.rds` + a `pi_models/` ensemble (+ a PI AUC/ROC eval) |
| **02 single-variant** | `marginal-scan` · `compute-ld-scores` · `calibrate` | the pheno bundle + the aGDS | `marginal_all.csv` + LD scores + a calibrated CSV |
| **03 SNV-set** | `prepare` (Stage 1) · `run-{gene,window,coding}` (Stage 2) | the pheno bundle + aGDS + `b_model` + `pi_models/` | per-unit flat tables (`results/…`) |
| **04 summary** | `aggregate` · `plots` | the Stage-2 per-unit tables | an aggregated table + top-hits + Manhattan/QQ |

`00` outputs are **cohort data assets** in a user `data_root` (a sibling-suffix
tree — see [Data lineage](#data-lineage-the-00-data-prep-tree)), not run `outputs/`.
`02`/`03`/`04` are **analysis** stages: they read those assets and write under the
run's `outputs/`.

## The three SNV-set styles

`region_type` (a config field) selects the analysis style; one **shared scan**
(`GLOWpipeline::run_scan_unit`) serves all three:

| region_type | unit | Stage-2 selector | per-unit output |
|---|---|---|---|
| `gene`   | gene                  | `--chr N`     | `results/glow_chr<N>.csv` |
| `window` | sliding window        | `--chunk K`   | `results/scan_chunk_<K>.csv` |
| `coding` | gene x coding category | `--chr N`     | `results-glow/glow_coding_chr<N>.csv` |

## Run organization

A **run** is a directory `runs/<name>/` holding a `config.R` and its sibling
`outputs/`. Every stage takes `--config <run>/config.R`; the scripts **derive**
`output_dir = <dir of --config>/outputs` and `run_name = <that dir's name>`.
`output_dir` is never declared in a config - that keeps config<->outputs
co-location structural and a run relocatable. `config.R` sources its
`base-config/*_base.R` and overrides only the cohort paths + what the run varies;
the run's diff is exactly its override block.

For a multi-stage cohort, the `runs/example/` layout below groups the per-stage
configs under one run directory with a shared cohort block:

```
runs/example/
  _cohort.R                     # the ONE place naming the (synthetic) data
  data-prep/config.R            # 00:   sources data_prep_base + _cohort
  training/config.R             # 01:   sources training_base + _cohort
  single-variant/config.R       # 02:   sources single_variant_base + _cohort
  snv-set/config.R              # 03/04 (gene):   sources snv_set_base + _cohort
  snv-set-window/config.R       # 03/04 (window): a region_type swap
  snv-set-coding/config.R       # 03/04 (coding): a region_type swap
```

Each per-stage `config.R` is `source(base); source(_cohort); <a few overrides>`;
`run_name`/`output_dir` derive from the per-stage dir (so `snv-set/` →
`snv-set/outputs/`). A bring-your-own cohort copies this layout and edits
`_cohort.R` (see [Bring your own cohort](#bring-your-own-cohort)).

## Data lineage (the 00-data-prep tree)

`00-data-prep` writes **cohort data assets** into a user `data_root` as a
**sibling-suffix tree** keyed by a `base_name` (resolved by
`00-data-prep/_dataprep_lib.R`):

```
<data_root>/<base_name>_plink/      raw PLINK (chr{chr}…)                [input]
<data_root>/<base_name>_gds/        plink-to-gds       → chr{chr}.gds
<data_root>/<base_name>_gds_favor/  annotate-favor     → <match>/{gds,csv}/   (the aGDS)
<data_root>/<base_name>_pcs/        compute-pcs        → pcs.rds / .csv
<data_root>/<base_name>_pheno/      assemble-pheno-covar → <run>_pheno_covar.rds
```

The analysis stages (`02`, `03`) point `gds_dir`/`gds_pattern` at the **aGDS** tree
(`<base_name>_gds_favor/<match>/gds/`) and `pheno_path` at the pheno bundle. Each
`00` step also drops a `README.md` + `provenance/` snapshot into its output dir, so
every derived dataset is self-documenting.

## End-to-end example (synthetic, out-of-the-box)

`data-example/` ships a small **synthetic** multi-chromosome cohort + a tiny bundled
FAVOR DB, all tracked, so the whole `00`→`04` chain runs on a fresh clone with **no
external data**. The outputs are *mechanically* valid but **statistically
meaningless** (it validates the template mechanics). The cohort is
shared-storage-independent (local `outputs/`, no `/project`).

```bash
cd /path/to/GLOWanalyses             # run from the GLOWanalyses directory
conda activate r_env                 # with GFisher + GLOWr + GLOWpipeline installed
bash run-example.sh                  # 02 → 03 → 04 over the pre-built assets, then the self-test
REBUILD=1 bash run-example.sh        # also regenerate the cohort + 00/01 intermediates
DRYRUN=1  bash run-example.sh        # print the command chain without running
```

`run-example.sh` runs each stage via its `--config` and finishes with
`data-example/check-example.R`, which asserts presence + schema + non-emptiness per
stage (and exits non-zero on any failure). What each stage emits for the example:

| Stage | Command (abbrev.) | Emits |
|---|---|---|
| 00 | `plink-to-gds` · `annotate-favor` · `compute-pcs` · `assemble-pheno-covar` | `example_gds/`, the scannable **aGDS** (16 PI + 3 coding nodes), `pcs.rds`, the pheno bundle |
| 01 | `b/{prepare,estimate}` · `pi/{prepare,annotate,train,evaluate}` | `b_model.rds`, a 5-model `pi_models/`, a PI AUC table + ROC |
| 02 | `marginal-scan` · `compute-ld-scores` · `calibrate` | `marginal_all.csv`, LD scores, `single_variant_all_calibrated.csv` |
| 03 | `prepare` · `run-gene --chr {21,22}` | `shared/01-*.rds`, `results/glow_chr{21,22}.csv` |
| 04 | `aggregate` · `plots` | `aggregated/glow_results_all.csv` + top-hits + Manhattan/QQ |

The committed cohort already includes the pre-built aGDS + B/PI models, so the
default `run-example.sh` (no `REBUILD`) skips `00`/`01` and runs only `02`→`04`.
**window** and **coding** are documented `region_type` swaps that run through the same scan —
`runs/example/snv-set-window/` and `runs/example/snv-set-coding/` smoke them
(`prepare` then one unit each). See `data-example/README.md` for the cohort shape,
the tracked-vs-generated split, and the one-command regeneration.

## Shared storage (optional)

GLOWanalyses templates default to **local** outputs (portable - they work without
a `/project` shared mount). To make a run's `outputs/` *born-shared* - physically
living off `/home` in a lab-readable shared root, with the repo keeping only a
symlink - **opt in per run** by setting `symlinked_shared_root` to a shared root
in that run's `config.R` (after it sources the base):

```r
symlinked_shared_root <- "/project/<lab>/<area>/GLOWpipeline-shared"
```

The stage then symlinks the run's `outputs/` into `<root>/<repo-relative-path>`
before populating it (honored by the `01-training` and `02-single-variant` run-org
templates and every `03-snv-set` template via `GLOWpipeline::prepare_scan_run`).
Leave the field `NULL` (the default) to keep `outputs/` a plain local directory.

`00-data-prep` outputs are *not* per-run born-shared - they are cohort data assets
under a user `data_root` (see [Data lineage](#data-lineage-the-00-data-prep-tree)),
shared at the `data_root` level if you place that root on a shared mount, not under
`runs/`. The bundled synthetic example leaves `symlinked_shared_root` unset, so it
runs against local outputs with no shared mount.

## End-to-end workflow (per run)

```bash
conda activate r_env
CFG=runs/<name>/config.R

# Stage 1 - shared prep (once): null model, snapshot, caches/chunks/STAAR nulls
Rscript 03-snv-set/prepare.R --config "$CFG"

# Stage 2 - scan (SLURM array; or run a unit directly to smoke-test)
#   gene:    array=1-22 over run-gene.sh   (per-chr)
#   window:  array=1-N  over run-window.sh (per-chunk; N from the chunk table)
#   coding:  array=1-22 over run-coding.sh (per-chr)
LOGS=runs/<name>/outputs/slurm-logs
sbatch --array=1-22 --output="$LOGS/%a_%j.log" --error="$LOGS/%a_%j.err" \
       slurm/run-gene.sh "$CFG"
# OR (recommended when the array exceeds the QOS MaxSubmitJobsPerUser cap):
#   drip-feed waves under CAP, auto-picking the array script + dimension from
#   the config's region_type; DRYRUN=1 previews without submitting.
CAP=90 bash slurm/submit-throttled.sh "$CFG"
# direct smoke (no SLURM):
Rscript 03-snv-set/run-gene.R --config "$CFG" --chr 22 --max-genes 20

# Stage 3 - aggregate + visualize
Rscript 04-summary/aggregate.R --config "$CFG"
Rscript 04-summary/plots.R     --config "$CFG"
```

All commands run from the **GLOWanalyses directory** (paths are
GLOWanalyses-root-relative). The Stage-2 entries take `--max-genes K`
(gene/coding) or `--max-windows K` (window) for fast smoke runs. `aggregate.R` takes `--top-k K` (gene) and `--allow-incomplete`
(window/coding partial-genome). `plots.R` takes `--drilldown` (gene SNV-level
per-top-hit CSVs; requires `write_evidence = TRUE`).

## Configuration reference (per stage)

Every stage sources its `base-config/*_base.R` (documented defaults, **no** cohort
paths) and the run config adds only what names real data + what the run varies. The
fields below are the ones a run typically sets; the base-config inline comments are
authoritative.

**`00` data-prep** (`data_prep_base.R`). *Required:* `data_root`, `base_name` (the
lineage anchor); `favor_db` (a FAVOR DB dir); the `assemble-pheno-covar` mapping —
`pheno_csv`, `pheno_id_col`, `outcome` (`list(col=|node=, map=)`), `covariates` (a
named list), `pcs` (`list(path=, id_col=, cols=)`). *Common:* `chroms`, the
`*_pattern` file patterns, `pc_source` (`gds`|`favor`), `n_pcs`, `pc_*` PCA knobs,
`favor_match_method`, `trait`, `sample_gds_chr`.

**`01` training** (`training_base.R`). *Required:* `b_raw_path` + `pi_raw_path` (the
known-SNP table); `favor_db`; `pi_control_source` (the FAVOR-annotated aGDS control
tree). *Common:* `b_method`/`b_qc_filters`/`b_column_mapping`; `pi_n_models`,
`pi_model_type`, `pi_controls_multiplier`, `pi_max_controls`, `pi_features`,
`pi_random_seed`, `pi_models_subdir`; `favor_features` (set to your DB's feature set
if it is not the full default).

**`02` single-variant** (`single_variant_base.R`). *Required:* `pheno_path` (the
`00` bundle), `gds_dir` + `gds_pattern` (the aGDS tree), `ld_scores_path` (a **shared**
genotype-only LD table). *Common:* `use_SPA`, `chunk_size`, `mac_cutoff`,
`missing_imputation`, `ld_window`/`ld_segment`, `calibration_method`.

**`03`/`04` SNV-set** (`snv_set_base.R`). *Required:* `gds_dir`, `gds_pattern` (with
the literal `{chr}` token), `pheno_path`, `pi_model_dir`, and **exactly one** B
source (`b_func` *or* `b_model_path`). *Common:* `region_type`, `chroms`,
`filter_spec` (`rare_maf_cutoff`/`variant_type`/`min_mac`/`min_variants`),
`ld_threshold`, `mac_threshold`, `collapse_method`, `pi_features`, `use_spa`
(NULL=auto), `calibration` (GC; default off), the `write_*` output-policy flags,
`output_format`, `primary_test`/`primary_test_glow`, `alpha`.
- *window-only:* `window_size`, `step_size`, `chunk_strategy`, `chunk_size_mb`,
  `merge_gap`, `cmac_cutoff`.
- *coding-only:* `categories` (a character vector of builtin names *or* a named list
  mapping each label to a builtin name / custom DNF clauses); the
  native-STAARpipeline comparison fields (`staar_modes`, `staar_bakein_glow`,
  `gene_num_in_array`, `Annotation_dir`, `staar_rv_num_cutoff`,
  `geno_missing_imputation`, `staar_genomewide_alpha`). `staar_enabled` toggles the
  STAAR comparator (a hard dependency only when TRUE); with it off, also set
  `staar_bakein_glow <- FALSE`.

`symlinked_shared_root` (all stages) is the shared-output opt-in (see below); leave
unset / `NULL` for local outputs.

## Bring your own cohort

1. **Prepare data assets (`00`).** Put per-chromosome PLINK (or GDS) under a
   `data_root`, point a data-prep config at it (`data_root`, `base_name`, `favor_db`,
   the pheno mapping), and run `plink-to-gds` → `annotate-favor` → `compute-pcs` →
   `assemble-pheno-covar`. (Already have GDS/aGDS? Skip `plink-to-gds`; set
   `gds_dir`/`gds_pattern` to your tree.)
2. **Train B + PI (`01`)** from a known-SNP table + your FAVOR-annotated control
   tree — or point `02`/`03` at a `b_model.rds` + `pi_models/` you already have.
3. **Analyze (`02`/`03`/`04`).** Copy the `runs/example/` layout, edit `_cohort.R`
   to name your data (`gds_dir`, `gds_pattern`, `pheno_path`, `pi_model_dir`, the B
   source, `pi_features`), and run each stage via its `--config`.

The `runs/example/` configs (synthetic) are the working template to copy: duplicate
the run directory under a new name, edit `_cohort.R` to point at your data, and run
each stage via its `--config`. Copying (rather than editing the shipped example in
place) keeps the bundled demo runnable and lets a run's `outputs/` derive from your
copy.


> Developed with AI assistance; see [AI-USE.md](AI-USE.md).
