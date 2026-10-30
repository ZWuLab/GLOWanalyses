# GLOWanalyses

Runnable, cohort-agnostic **workflow templates** for the GLOW methodology. This
is a *scripts directory, not an R package* (no `DESCRIPTION`/`NAMESPACE`): copy a
run, point it at your data via a config, and run. The statistical methods live in
`GLOWr`; the orchestration + summary live in `GLOWpipeline`; these scripts are
the thin, numbered glue that calls them.

```
00-data-prep/  PLINK→GDS · FAVOR annotation (→aGDS) · PCs · pheno/covar assembly
01-training/   B (effect-size dist.) + PI (variant-importance score) training
02-single-variant/  per-variant scan · in-sample LD scores · GC calibration
03-snv-set/    Stage 1 (prepare.R) + thin Stage-2 entries (run-gene/window/coding/noncoding.R)
               + optional native-STAARpipeline comparison engines (run-staar-coding/window.R)
04-summary/    Stage 3: aggregate.R (region_type-aware) + plots.R (Manhattan/QQ/lambda/drilldown)
base-config/   {data_prep,training,single_variant,snv_set}_base.R - documented defaults (no cohort paths)
data-example/  a self-contained synthetic cohort + generator + self-test (runs out-of-the-box)
runs/<name>/   per-run config.R (sources the base, names the data) + derived outputs/
slurm/         array-job templates: 00 (run-plink-to-gds/run-annotate-favor.sh), 02 (run-marginal.sh),
               03 (run-gene/window/coding/noncoding.sh) + submit-throttled.sh; _job_lib.sh activates R inside each job
cloud/         the official Docker image recipes (docker/) and the Terra workflow glow_scan.wdl (terra/):
               the 02-04 stages on Terra or any WDL engine, over these same scripts (see cloud/README.md)
```

All five stages (`00`→`04`) are thin `--config` wrappers over the packaged
`GLOWr` / `GLOWpipeline` methods. The repository ships a **runnable,
self-contained synthetic example** (`data-example/` + `runs/example/`) — run
`bash run-example.sh` from the GLOWanalyses directory and the whole chain
executes on a tiny tracked cohort (see the
[End-to-end example](#end-to-end-example-synthetic-out-of-the-box) walkthrough
below). The four SNV-set types (gene, sliding window, coding masks and the gene-centric
non-coding masks) run through the same scan; see "How GLOW forms SNV sets".

## Install + run convention

These scripts call the three GLOW packages at runtime; **install them first**
(they are separate repos), then run every command **from the GLOWanalyses
directory** (paths are GLOWanalyses-root-relative):

```r
# in R, once: install the runtime dependencies (GFisher -> GLOWr -> GLOWpipeline)
# e.g. remotes::install_github("ZWuLab/GLOWr") etc., or from local source trees
```

Beyond the three GLOW packages (and what they import), the templates attach these
Bioconductor and CRAN packages themselves, so install them too:

| Package | Needed by |
| --- | --- |
| `SeqArray`, `gdsfmt` | every stage that reads a GDS or aGDS |
| `SNPRelate` | `00-data-prep/compute-pcs.R`: the LD pruning and the PCA that `GLOWr::compute_pcs_gds()` delegate to it |
| `arrow` | `00-data-prep/annotate-favor.R` with a FAVOR 2.0 Parquet database |
| `STAAR` | the STAAR comparator inside the GLOW scan (`staar_enabled <- TRUE`, any region type) |
| `STAARpipeline`, `SeqVarTools` | the native STAARpipeline runners `03-snv-set/run-staar-*.R` only |

```bash
cd /path/to/GLOWanalyses          # all stage/config paths resolve from here
conda activate r_env              # or any R with the 3 packages installed
```

**Inside SLURM jobs** the `slurm/run-*.sh` templates activate the R environment
themselves (a batch job inherits only environment variables, not an active
shell environment): set `GLOW_CONDA_ENV` to your conda env (default `r_env`),
or `GLOW_RSCRIPT` to an `Rscript` path (e.g. from `module load R`), before
`sbatch` — e.g. `export GLOW_CONDA_ENV=GLOW`. See `slurm/_job_lib.sh`.

## Getting and updating your copy

Clone the repository rather than downloading a zip: a clone updates in place and
always tells you which version you are running.

```bash
git clone https://github.com/ZWuLab/GLOWanalyses.git
cd GLOWanalyses
cat .release-source        # the version + the source commit this copy was built from
```

To keep updates painless:

- **Keep your own work in `runs/<name>/`** — your configs and their `outputs/`
  (which the shipped `.gitignore` ignores) — and your data outside the clone (the
  `data_root`). `git status` then shows only your files as untracked, and
  `git pull` never touches them.
- **Do not edit the shipped scripts in place.** Copy one under a new name (e.g.
  `slurm/my-run-gene.sh`) and edit the copy, so an update cannot collide with it;
  cluster-specific settings already have knobs (`GLOW_CONDA_ENV`, `GLOW_RSCRIPT`,
  `sbatch` flags on the command line).
- **Update with `git pull`**, or pin a version with `git checkout v0.1.1`. If git
  refuses because one of your untracked files has the same path as a newly shipped
  file, move your file aside and pull again.
- **Never run `git clean` here** — it deletes untracked files, i.e. your configs and
  logs. `git pull`, `git checkout` and `git status` are safe.

Do not commit or push anything to the repository; report problems or
suggestions through the repository's issues.

## The five stages (`00`→`04`)

Each stage is run from the GLOWanalyses directory as `Rscript
<stage>/<script>.R --config <run>/config.R [unit selector]`. The
producer/consumer chain:

| Stage | Scripts | Reads | Writes |
|---|---|---|---|
| **00 data-prep** | `plink-to-gds` · `annotate-favor` · `compute-pcs` · `assemble-pheno-covar` | raw PLINK + a FAVOR DB + a pheno table | the cohort **data-lineage tree** (`<base>_gds`, `_gds_favor` (aGDS), `_pcs`, `_pheno`) |
| **01 training** | `b/{prepare,estimate}` · `pi/{prepare,annotate,train,evaluate}` | a known-SNP table + the FAVOR DB + the aGDS control tree | `b_model.rds` + a `pi_models/` ensemble (+ a PI AUC/ROC eval) |
| **02 single-variant** | `marginal-scan` · `compute-ld-scores` · `calibrate` | the pheno bundle + the aGDS | `marginal_all.csv` + LD scores + a calibrated CSV |
| **03 SNV-set** | `prepare` (Stage 1) · `run-{gene,window,coding,noncoding}` (Stage 2) | the pheno bundle + aGDS + `b_model` + `pi_models/` | per-unit flat tables (`results/…`) |
| **04 summary** | `aggregate` · `plots` | the Stage-2 per-unit tables | an aggregated table + top-hits + Manhattan/QQ |

`00` outputs are **cohort data assets** in a user `data_root` (a sibling-suffix
tree — see [Data lineage](#data-lineage-the-00-data-prep-tree)), not run `outputs/`.
`02`/`03`/`04` are **analysis** stages: they read those assets and write under the
run's `outputs/`.

## How GLOW forms SNV sets

`region_type` (a config field) selects how the SNV sets are formed, and one **shared scan**
(`GLOWpipeline::run_scan_unit`) serves all four types. A scan unit is one chromosome (for
the `gene`, `coding` and `noncoding` types) or one genome chunk (for `window`). A unit holds
many regions, and GLOW tests one SNV set per region. The set is formed in two steps that the
code keeps apart. The **region** says which variants are candidates. The **filter**
(`filter_spec`: the QC label, the variant type, the MAF and MAC cutoffs, the minimum set size)
says which candidates stay. The four types differ in how the region is defined:

| `region_type` | one set per | candidates | how a variant is tied to the set | Stage-2 selector | per-unit output |
|---|---|---|---|---|---|
| `gene` | gene, or segment of a gene above `gene_max_records` | every record in the gene's or segment's span | by position | `--chr N` | `results/glow_chr<N>.csv` |
| `window` | sliding window | every record in the window | by position | `--chunk K` | `results/scan_chunk_<K>.csv` |
| `coding` | gene and coding mask | the records in the gene's span whose coding annotation matches the mask | by position, then by annotation | `--chr N` | `results-glow/glow_coding_chr<N>.csv` |
| `noncoding` | gene and non-coding mask | the records anywhere on the chromosome whose annotation names the gene and matches the mask | by annotation only | `--chr N` | `results-glow/glow_noncoding_chr<N>.csv` |

**`gene`.** One set per gene of the bundled protein-coding gene table (18,445 genes), or of
your own list (`resolved_gene_list`). On whole-genome sequencing data a gene with more records than
`gene_max_records` is tiled into segments of `gene_segment_bp`, each tested as its own row
labeled `<gene>_seg<k>` (see "Cores and large genes").

**`window`.** One set per sliding window of `window_size` base pairs every `step_size`,
over the whole chromosome, so non-coding regions are covered by position without any mask.

**`coding`.** One set per gene and coding mask, the records in the gene's span whose
coding annotation matches the mask. The built-in masks follow STAARpipeline: `plof`,
`plof_ds`, `missense`, `disruptive_missense`, `synonymous`, `ptv`, `ptv_ds` (`categories`).

**`noncoding`: the gene-centric non-coding masks.** One set per gene and mask, the records
that the annotation links to the gene wherever they lie on the chromosome, as in
STAARpipeline's gene-centric non-coding analysis. A variant upstream of a gene lies outside
the gene body, so gene coordinates cannot locate it. The link comes from the FAVOR annotation
of the record. The eight masks (`noncoding_categories`):

| mask | gene table | a record is selected when | the gene it is linked to |
|---|---|---|---|
| `upstream` | protein-coding | its GENCODE category has the token `upstream` | every gene the upstream group of its GENCODE info names |
| `downstream` | protein-coding | its GENCODE category has the token `downstream` | every gene the downstream group names |
| `UTR` | protein-coding | its GENCODE category has `UTR3` or `UTR5` | every gene its UTR groups name |
| `promoter_CAGE` | protein-coding | its `cage_tc` field is non-empty and it lies within 3 kb of a TSS | the first gene its GENCODE info names |
| `promoter_DHS` | protein-coding | its `rdhs` field is non-empty and it lies within 3 kb of a TSS | the first gene its GENCODE info names |
| `enhancer_CAGE` | protein-coding | its `cage_tc` and `genehancer` fields are non-empty | the first connected gene of its GeneHancer element |
| `enhancer_DHS` | protein-coding | its `rdhs` and `genehancer` fields are non-empty | the first connected gene of its GeneHancer element |
| `ncRNA` | the bundled ncRNA gene table (21,104 names) | its GENCODE category has `ncRNA_exonic` or `ncRNA_splicing` | every gene its GENCODE info names |

A record can therefore join several genes' sets under `upstream`, `downstream`, `UTR` and
`ncRNA`. Under the promoter and enhancer masks it joins one gene's set. The terms:

- GENCODE: the gene annotation whose category (upstream, downstream, UTR, exonic, intronic
  and so on) and gene names FAVOR carries per record.
- UTR: an untranslated region of a transcript. ncRNA: a non-coding RNA gene.
- CAGE and DHS: experimental marks of active promoters and of open chromatin, recorded by
  FAVOR as the `cage_tc` and `rdhs` fields.
- GeneHancer: a database of enhancer elements with the genes they are connected to,
  recorded as the `genehancer` field.
- TSS: a transcription start site.

The five aGDS fields the masks read are `genecode_comprehensive_category`,
`genecode_comprehensive_info`, `genehancer`, `cage_tc` and `rdhs`, which
`00-data-prep/annotate-favor.R` writes by default. Both FAVOR layouts are read, the FAVOR v1
strings and the FAVOR 2.0 serialization of the annotator. The FAVOR 2.0 layout does not keep
the group boundary of a record with two categories, and such a record is then linked to
every gene it names under each of its masks. The two promoter masks take
the transcription start sites from `TxDb.Hsapiens.UCSC.hg38.knownGene` (with
`GenomicFeatures`), genome-wide, so a variant can lie in gene A's window and be linked to
gene B by its annotation. Your aGDS must be on hg38 for these two masks. We call the eight
masks the gene-centric non-coding masks, following STAARpipeline. Four of them (`upstream`,
`downstream`, `UTR`, `ncRNA`) select by GENCODE category. The promoter and enhancer masks
select by the CAGE, DHS and GeneHancer annotations instead, so they can include exonic
variants. The run's `filter_spec` applies to every mask without an added annotation clause,
and the sets are LD-pruned, collapsed and tested exactly as the other types.

GLOW reads the annotation by its meaning rather than by STAARpipeline's string operations,
and a set can differ from STAARpipeline's in four cases. A record annotated
`upstream;downstream` joins both masks. A UTR record naming two genes joins both genes'
sets. A hyphenated gene name is kept whole in the promoter masks. An ncRNA record with
several gene groups links them all. We validated this on chromosome 22 of the 1000 Genomes
Project (3,202 samples, 427 protein-coding and 483 ncRNA genes). The sets formed with
STAARpipeline 0.9.8's own string operations (`rule = "staar"` of
`GLOWr::define_regions_noncoding()`) matched the native selection variant for variant, and
the default rule differed only in those four cases. The per-unit table has the coding
table's shape, with `category` the mask. A row's `start`, `end` and `midpoint` describe the
span of the candidates, which is wide for an enhancer far from its gene. The aggregation
orders a gene's masks together and computes the Bonferroni cutoff per mask and across the
configured masks, ncRNA included when it is selected, from the tested cells
(`aggregated/noncoding_thresholds.csv`, every cutoff with its denominator). Beside them it
prints, for reference, the constants of the STAARpipeline paper: 0.05 / (20,000 x 7) for the
seven masks of protein-coding genes, and 0.05 / 20,000 for the ncRNA genes (Li, Z., Li, X.,
et al. 2022, Nature Methods 19,
1599-1611, doi 10.1038/s41592-022-01640-x). In prose, "non-coding" names a class of variants
or regions, which the sliding-window scan covers by position. "Gene-centric non-coding
masks" names this family.

## Run organization

A **run** is a directory `runs/<name>/` holding a `config.R` and its sibling
`outputs/`. Every stage takes `--config <run>/config.R`; the scripts **derive**
`output_dir = <dir of --config>/outputs` and `run_name = <that dir's name>`.
`output_dir` is never declared in a config - that keeps config<->outputs
co-location structural and a run relocatable. `config.R` sources its
`base-config/*_base.R` and overrides only the cohort paths + what the run varies;
the run's diff is exactly its override block.

**Three layers of configuration, and which ones you edit.** A stage's `config.R` is
`source(base); source(_cohort); <overrides>`, so a field set later wins:

1. `base-config/*_base.R`: the shipped defaults, with the documentation of every field in
   their inline comments. **Do not edit these files.** They belong to the repository, and
   `git pull` overwrites them with the next release.

2. `runs/<name>/_cohort.R`: the one place that names your data (paths, file-name patterns,
   the feature set). Yours to write.

3. `runs/<name>/<stage>/config.R`: what this run varies (the region type, the variant
   filter, the STAAR comparator, calibration). Yours to write.

Put every change of yours in layers 2 and 3. A field set there replaces the same field of
`base-config/`, because it is sourced afterwards.

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

**Where outputs go.** Two places, by stage:

| Stage | Writes to | What |
| --- | --- | --- |
| `00` data-prep | the `data_root` tree that `_cohort.R` names (next section) | the cohort data assets: GDS, aGDS, PCs and the phenotype bundle, each with its `provenance/` folder |
| `01` to `04` | `runs/<name>/<stage>/outputs/`, created next to the stage's `config.R` | the models, scans, results, logs and plots of that run |

The analysis stages never write under `data_root`, and `00` writes nothing under `runs/`.

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

Two optional fields move the annotation step off this layout: `favor_input_gds_dir`
reads another GDS tree than `<base_name>_gds/`, and `favor_output_dir` writes one flat
tree `<dir>/{gds,csv,provenance,logs}` (no `<match>` level), which `compute-pcs` and
`assemble-pheno-covar` then read. Use them to place a re-annotation beside, not over,
an earlier tree. Left unset, the layout above applies unchanged.

The analysis stages (`02`, `03`) point `gds_dir`/`gds_pattern` at the **aGDS** tree
(`<base_name>_gds_favor/<match>/gds/`) and `pheno_path` at the pheno bundle. Each
`00` step also drops a `README.md` + `provenance/` snapshot into its output dir, so
every derived dataset is self-documenting.

At whole-genome scale run the two per-chromosome `00` steps as **SLURM arrays**
(one task per chromosome; FAVOR annotation is the heaviest step — ~32 GB peak and
~15 min per chromosome measured against the FAVOR Essential DB), then the two
serial steps:

```bash
CFG=runs/<name>/data-prep/config.R
LOGS=runs/<name>/data-prep/outputs/slurm-logs; mkdir -p "$LOGS"
J1=$(sbatch --parsable --array=1-22 \
            --output="$LOGS/gds_chr%a_%j.log" --error="$LOGS/gds_chr%a_%j.err" \
            slurm/run-plink-to-gds.sh "$CFG")
sbatch --dependency=afterok:$J1 --array=1-22 \
       --output="$LOGS/favor_chr%a_%j.log" --error="$LOGS/favor_chr%a_%j.err" \
       slurm/run-annotate-favor.sh "$CFG"
# after the FAVOR array finishes (serial, no array):
Rscript 00-data-prep/compute-pcs.R          --config "$CFG"
Rscript 00-data-prep/assemble-pheno-covar.R --config "$CFG"
```

Each array task writes its own `provenance/config_snapshot_chr<N>.rds`; `#SBATCH`
resources are dev-cluster defaults — override them on the `sbatch` line.

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

## Running on Terra (WDL + Docker)

The scan stages also run in the cloud, without installing anything. `cloud/terra/glow_scan.wdl`
is a WDL workflow whose tasks run the templates of `02-single-variant`, `03-snv-set` and
`04-summary` inside the official image `ghcr.io/zwulab/glowanalyses:<version>`. The image
carries R, GFisher, GLOWr, GLOWpipeline and this tree at matching versions.
`region_type = "single_variant"` gives the per-variant scan with its calibration, and `gene`,
`window`, `coding` or `noncoding` give the SNV-set scan with its aggregation, calibrated by the
single-variant run's calibration report (one file) or its two files. The analysis settings
carry the names of the configuration fields below (the file inputs are translated into the
template paths, see `cloud/README.md`), every output is one named table (a Terra data-table
column each, with ready-made configuration files per run type in `cloud/terra/configs/`), and
on the bundled example the tables equal this route's. Import the workflow into a Terra
workspace, or run it on any machine
with [miniwdl](https://miniwdl.readthedocs.io) and Docker or Apptainer.
`bash cloud/terra/run-example-wdl.sh` runs the synthetic example through it and compares with
`run-example.sh`. Details, sizing, the defaults to check and the image's contents:
`cloud/README.md`.

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
#   noncoding: array=1-22 over run-noncoding.sh (per-chr; the eight masks per gene)
LOGS=runs/<name>/outputs/slurm-logs
sbatch --array=1-22 --output="$LOGS/%a_%j.log" --error="$LOGS/%a_%j.err" \
       slurm/run-gene.sh "$CFG"     # coding: slurm/run-coding.sh; noncoding: slurm/run-noncoding.sh
# OR (recommended when the array exceeds the QOS MaxSubmitJobsPerUser cap):
#   drip-feed waves under CAP, auto-picking the array script + dimension from
#   the config's region_type; DRYRUN=1 previews without submitting.
CAP=90 bash slurm/submit-throttled.sh "$CFG"
# direct smoke (no SLURM):
Rscript 03-snv-set/run-gene.R --config "$CFG" --chr 22 --max-genes 20
# the other types take the same arguments with their own entry script, e.g. for
# region_type <- "noncoding" (prepare.R once, then one chromosome at a time):
Rscript 03-snv-set/run-noncoding.R --config "$CFG" --chr 22 --max-genes 20   # caps every mask's table

# OPTIONAL - the STAAR comparator inside the GLOW scan (gene, window, coding or noncoding; needs
# the STAAR package): set `staar_enabled <- TRUE` in the config BEFORE Stage 1.
# prepare.R then fits the STAAR null model on GLOW's phenotype and covariates, and
# every gene / window / coding / noncoding cell also carries STAAR::STAAR() on GLOW's own variant
# set: STAAR_O, ACAT_O and the SKAT / Burden / ACAT-V family omnibi at both beta
# weights, in the same result tables as the GLOW tests. The comparator hands STAAR the
# genotype matrix as the GDS stores it (the reference-allele dosage) and STAAR flips it to
# the minor allele inside; STAARpipeline flips before its call, so on an identical variant
# set the two STAAR-O values agree to about 1e-8 relative, not to machine precision.

# OPTIONAL - native STAARpipeline comparison (needs the STAARpipeline package;
# submit once per SPA mode; window mode needs staar_native = TRUE in the config
# before Stage 1 so prepare.R fits the native null models):
#   coding: array over run-staar-coding.sh (per gene-group; N in the .sh header)
#   window: array=1-N over run-staar-window.sh (per-chunk, the GLOW grid exactly)
sbatch --array=1-"$N" --output="$LOGS/sw_off_%a_%j.log" --error="$LOGS/sw_off_%a_%j.err" \
       slurm/run-staar-window.sh "$CFG" off

# Stage 3 - aggregate + visualize
Rscript 04-summary/aggregate.R --config "$CFG"
Rscript 04-summary/plots.R     --config "$CFG"
```

All commands run from the **GLOWanalyses directory** (paths are
GLOWanalyses-root-relative). The Stage-2 entries take `--max-genes K`
(gene/coding/noncoding; for noncoding it caps every mask's table) or `--max-windows K` (window) for fast smoke runs. `aggregate.R` takes `--top-k K` (gene) and `--allow-incomplete`
(window/coding/noncoding partial-genome). `plots.R` takes `--drilldown` (gene SNV-level
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
`favor_match_method`, `favor_db_format`, `favor_release`, `favor_rsid_policy`,
`favor_input_gds_dir`, `favor_output_dir`, `trait`, `sample_gds_chr`.

*FAVOR annotation.* `favor_db` is a FAVOR v1 CSV directory (the `chr*_*.csv` chunks) or a FAVOR
2.0 directory (one `chromosome_<N>.parquet` per chromosome). You obtain the FAVOR 2.0 files from
the FAVOR team. The field `favor_db_format` is `"auto"` by default and detects which one it is,
and the Parquet format needs the R `arrow` package. The setting `favor_match_method = "exact"`
looks each variant up by its own CHR-POS-REF-ALT. The setting `"flexible"` also tries, for SNVs,
the key with REF and ALT swapped and the key on the opposite strand, which genotyping-chip data
converted through PLINK can need. Your genotypes and alleles are never changed. Only the lookup
key is. Every variant gets an outcome,
`match_tier` in the CSV and `annotation/info/favor_match_tier` in the aGDS: `exact`, `swapped`,
`flipped`, `flipped_swapped`, `position_only` (a key without alleles), or a reason it was not
annotated (`rsid_conflict`, `unmatched_alt`, `unmatched_ref`, `uncovered`, `no_database`,
`unsupported_allele`). A swapped or flipped outcome says which lookup found the row. It is not a
check of the strand. The evidence for a match is the rsID: FAVOR carries dbSNP's rsID on the
specific REF-ALT row dbSNP knows, so every matched variant also records whether the chip's rsID
(the GDS `annotation/id`) equals the FAVOR row's, in `rsid_check` (`same`, `differs`,
`chip_none`, `favor_none`) and `annotation/info/favor_rsid_check`. The field `favor_rsid_policy`
is `"require"` by default: a swapped or flipped match that the rsID does not confirm is withheld
as `rsid_conflict`, so a rewritten key is accepted only when dbSNP confirms it or the data
carries no rsID to check. `"record"` keeps every match and only records the check, for
diagnostics. Position-only keys are annotated with the average over the SNV rows at
the position. The field `favor_release` is a label recorded in the aGDS, since the FAVOR files
do not state their release. The vignette `vignette("favor-annotation", package = "GLOWr")`
explains the matching rule, the tiers and the rsID evidence.

**`01` training** (`training_base.R`). *Required:* `b_raw_path` + `pi_raw_path` (the
known-SNP table); `favor_db`; `pi_control_source` (the FAVOR-annotated aGDS control
tree). *Common:* `b_method`/`b_qc_filters`/`b_column_mapping`; `pi_n_models`,
`pi_model_type`, `pi_controls_multiplier`, `pi_max_controls`, `pi_features`,
`pi_random_seed`, `pi_models_subdir`; `favor_features` (set to your DB's feature set
if it is not the full default).

**`02` single-variant** (`single_variant_base.R`). *Required:* `pheno_path` (the
`00` bundle), `gds_dir` + `gds_pattern` (the aGDS tree), `ld_scores_path` (a **shared**
genotype-only LD table). *Common:* `use_SPA`, `chunk_size`, `mac_cutoff`,
`missing_imputation`, `n_cores` (forked workers per chromosome scan; default 1),
`ld_window`/`ld_segment`, `calibration_method`.

**`03`/`04` SNV-set** (`snv_set_base.R`). *Required:* `gds_dir`, `gds_pattern` (with
the literal `{chr}` token), `pheno_path`, `pi_model_dir`, and **exactly one** B
source (`b_func` *or* `b_model_path`). *Common:* `region_type`, `chroms`,
`filter_spec` (`rare_maf_cutoff`/`variant_type`/`min_mac`/`min_variants`; the shipped `0.5`,
`"SNV"`, `1L`, `1L` keep the whole frequency spectrum and one-variant sets, GLOW's design rather
than STAARpipeline's rare-variant conventions; `cloud/README.md`, "Defaults to check before a
real run", has the override lines for either route),
`ld_threshold`, `mac_threshold`, `collapse_method`, `n_cores` (forked workers per scan
unit; default 1, see "Cores and large genes" below; a Stage-2 script's `--n-cores N` overrides
the snapshot's value for that unit), `pi_features`, `use_spa`
(NULL=auto), `debug_errors` (FALSE; TRUE puts the failing call and the stack of function
names into the `message` of an error row), `calibration` (GC; default off: the two files of a single-variant run, or its
`report`, the one-file `single_variant_calibration.tsv` that `02-single-variant/calibrate.R`
writes; Stage 1 refuses a report of another phenotype bundle), the `write_*` output-policy flags,
`output_format` (`csv`; `fst` is accepted for window runs only, because the gene, coding and
noncoding tables are read as CSV, and Stage 1 refuses it elsewhere), `top_k` (rows of
`glow_top_hits.csv` in gene runs; default 100), `primary_test`/`primary_test_glow`, `alpha`, `staar_enabled` (the STAAR
comparator on GLOW's own variant sets, for gene, window, coding and noncoding runs alike: Stage 1 fits
the STAAR null model on GLOW's phenotype and covariates, and every result row also carries
`STAAR_O`, `ACAT_O` and the six family omnibi, named `STAAR_O_glowG`/`ACAT_O_glowG` in the
coding aggregates; STAAR is a hard dependency only when TRUE), `staar_anno_features` (its
annotation-weight columns; `NULL` = the PI features), `write_staar_detail`.
- *gene-only:* `gene_max_records` (default 10,000) and `gene_segment_bp` (default 10,000):
  a gene with more GDS records in its span than `gene_max_records` is tiled into
  non-overlapping segments of `gene_segment_bp` from the gene start, each reported as its
  own row labeled `<gene>_seg<k>` (see "Cores and large genes" below).
- *window-only:* `window_size`, `step_size`, `chunk_strategy`, `chunk_size_mb`,
  `merge_gap`, `cmac_cutoff`.
- *coding-only:* `categories` (a character vector of builtin names *or* a named list
  mapping each label to a builtin name / custom DNF clauses); the
  native-STAARpipeline comparison fields (`staar_modes`, `staar_bakein_glow`,
  `gene_num_in_array`, `Annotation_dir`, `staar_rv_num_cutoff`,
  `geno_missing_imputation`, `staar_genomewide_alpha`). With `staar_enabled` off, Stage 1
  records `staar_modes` as empty and `staar_bakein_glow` as off, so a GLOW-only coding run
  needs nothing else.
- *noncoding-only:* `noncoding_categories` (the masks to run; default all eight, see "How
  GLOW forms SNV sets"). The per-chromosome table is
  `results-glow/glow_noncoding_chr<N>.csv`, and `aggregate.R` writes
  `glow_noncoding_all.csv`, `noncoding_merged.csv`, the two significance tables,
  `noncoding_thresholds.csv` (every cutoff with its denominator, STAARpipeline's
  constants beside it) and `noncoding_completeness.csv`. The two promoter masks need
  the `TxDb.Hsapiens.UCSC.hg38.knownGene` and `GenomicFeatures` packages.
  The per-chromosome table's columns: `label` (the gene), `chr`, `category` (the mask);
  `start`, `end` and `midpoint` (the first and last candidate position and their midpoint,
  the candidate span, which is wide for an enhancer far from its gene); the count chain of
  the set as in every GLOW table (`n_variants`, the candidates in the file;
  `n_after_annotation`, after the QC and type settings; `n_after_filter`, after the MAF and
  MAC cutoffs; `n_after_ld` and `n_after_collapse`, after the LD step and the collapse;
  `cMAC`; `is_single`); the GLOW tests, `GLOW_Omni` the primary one; and with
  `staar_enabled` the comparator's `STAAR_O`, `ACAT_O` and six family omnibi (renamed
  `STAAR_O_glowG` and `ACAT_O_glowG` in the aggregates). `noncoding_thresholds.csv` has one
  row per method and mask and one per method for `all`, the configured masks together:
  `n_tested`, the computed `cutoff`, its `cutoff_rule` (`alpha / n_tested` per mask,
  `family_alpha / n_all` across) and the reference constant with its derivation.
  `noncoding_completeness.csv` has one row per requested chromosome: its `status` (present,
  empty results file, no results file), the cells and genes, the tested and missing primary
  p-values, and the cells attempted, ok, skipped and errored from the per-mask summary the
  scan writes beside each table (`glow_noncoding_chr<N>_summary.csv`; the errors beside it in
  `glow_noncoding_chr<N>_errors.csv`, zero rows on a clean run).

`symlinked_shared_root` (all stages) is the shared-output opt-in (see below); leave
unset / `NULL` for local outputs.

### Cores and large genes

The SLURM templates run one single-core task per chromosome or chunk, which is the right
layout on a cluster with a shared filesystem. On a cloud platform such as Terra, where a task
is a virtual machine that first copies its inputs, the economical layout is one task per
chromosome with the machine's cores filled from inside R. Set `n_cores` in the run config to
that number: `run_scan_unit()` then splits the chromosome's regions into blocks and runs them
in forked workers, each with its own GDS handle, and the single-variant scan does the same with
its variant chunks. The written tables are identical to the single-core ones. Size the task's
memory as `n_cores` times the per-region peak, and set `OPENBLAS_NUM_THREADS=1` on a machine
with a threaded BLAS so the workers do not oversubscribe it. For scale: the `noncoding` unit of
chromosome 22 of the 1000 Genomes Project (3,202 samples, eight masks, the STAAR comparator on,
four workers) had a peak R heap of 1.3 GB in the parent process and a job-wide peak of 8.4 GB;
on a 16 GB node start with two workers and raise `n_cores` once a large chromosome has run. On a SLURM cluster, raise
`--cpus-per-task` and `--mem` together with `n_cores`. Forking is not available on Windows.

The regions of a chromosome are not handed to the workers in genomic order. Genes with many
records sit next to each other, and a block of neighbors would keep one worker busy after the
others have finished. Instead, `run_scan_unit()` reads every region's record count from the
file's index and takes the squared count as the region's cost. A region whose cost reaches the
chromosome's total cost divided by 8 times `n_cores` runs in a block of its own. The remaining
regions are packed, in decreasing cost, into blocks of at most 25 regions whose summed cost
stays within that threshold, and the blocks are dispatched heaviest first. The aim is a wall
time close to the total work divided by the cores. For the gene, window and coding types the
record count is read from the file's index over the region's span; for the non-coding masks
it is the region's candidate count. The budget and the blocks are formed per regions table:
once per unit for a gene or window run, once per chromosome for the coding categories, which
share the gene table, and once per mask for the non-coding type. On chromosome 22 of the 1000
Genomes Project (3,202 samples) with twenty workers on a node allocated to the job alone, this
order took 55 and 56 percent of the time of the genomic-order blocks in two repeats; with four
workers the two orders were within about 1 percent of each other. The
per-chromosome log shows one scheduler line per regions table (regions, blocks, how many run
alone, the heaviest region) and one line per block with its region count, its heaviest record
count and its wall time, written after all of the table's blocks have returned, in dispatch
order, so you can see where the time went. A gene list with a repeated label stops the run,
at any `n_cores`, because the per-region output files are keyed by label.

On whole-genome-sequencing data a gene body can hold tens of thousands of variant records, and
the per-gene cost and memory grow with the square of the variant count. A gene run therefore
tiles every gene with more GDS records in its span than `gene_max_records` into non-overlapping
segments of `gene_segment_bp` base pairs, counted from the file's positions before any filter,
and tests and reports each segment as its own row labeled `<gene>_seg<k>`, with no gene-level
combination. This is how the GLOW methodology paper handled genes above 1,000 variants, with
10 kb windows. The defaults, 10,000 records and 10 kb, leave genotyping-chip data untouched.
Lower `gene_max_records` for a large cohort (the genotype matrix of one gene is 8 bytes times
the sample count times the variant count), or set it to `Inf` to never segment.

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
