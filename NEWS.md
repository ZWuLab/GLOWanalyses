# GLOWanalyses 0.3.1

Requires GLOWr 0.3.3 and GLOWpipeline 0.2.1 (unchanged; the image's packages are those of
0.3.0). This release revises the Terra workflow on the feedback of Alberto Brusati (Landers
lab, UMass Chan Medical School), who runs it on Terra, and extends the templates behind it.

- **Named outputs of `cloud/terra/glow_scan.wdl`.** Every table is its own workflow output
  (`combined_results`, `top_hits`, `sig_windows`, `loci`, the coding and noncoding tables,
  `lambda_by_category`, `plots`, `run_info`, `per_chrom_results`, `per_chrom_errors`,
  `per_chrom_summaries`, `per_chrom_logs`, the `*_resources` files; the single-variant tables
  and `calibration_report_out`), so that a Terra data table holds one column per table. The
  earlier `scan_results`, `scan_logs`, `aggregated`, `marginal_all` and `single_variant_logs`
  outputs are renamed or split accordingly; `aggregated_tar` and `shared_artifacts` stay.
- **Terra configuration files, `cloud/terra/configs/`.** One inputs file and one outputs file
  per run type (`GLOWScan_<type>.inputs.json`, `GLOWScan_<type>.outputs.json`), in the format
  Terra's own "Download json" writes, binding the inputs to the data-table columns `gds_files`,
  `chromosomes`, `pheno_bundle`, `pi_models_tar`, `b_model_rds` and `calibration_report`, and
  writing each output to a `glow_<type>_<output>` column. The single-variant run writes its
  report to `calibration_report`, which the set-based runs read. `cloud/README.md` has the
  column table.
- **The workers come from the machine.** Each scan task uses the vCPUs it runs on (`nproc`) as
  `n_cores`, not the `cpu` request, because a memory request can make the platform give more
  vCPUs than asked and those cores are paid for; the new `n_cores` input overrides it. The four
  Stage-2 scripts gain `--n-cores N`, which overrides the Stage-1 snapshot's value for that
  unit on either route.
- **A resource monitor in every task, `cloud/terra/task-monitor.sh`.** Every 10 minutes and at
  exit, the memory in use and its peak, the load average and the vCPUs, to the task's
  `resources.txt` and its log, the exit line also on failure. It is how `mem_gb_scan` is sized
  per task kind.
- **Calibration from one file.** `02-single-variant/calibrate.R` also writes
  `single_variant_calibration.tsv`, a key/value report with the factor at full precision, its
  method, the regression's size, the lambdas, the phenotype bundle's checksum and the versions.
  `03-snv-set/prepare.R` takes it as `calibration$report` (the workflow input
  `calibration_report`) in place of `single_variant_scan` and `ld_scores`, and refuses a
  report whose `pheno_bundle_md5` is not the run's bundle. A run calibrated from the report
  equals one calibrated from the two files, byte for byte.
- **Stage-1 checks for a non-coding run.** `prepare.R` checks the five annotation fields the
  masks read on every chromosome, and the two packages the promoter masks need, before any
  scan starts, so a wrong aGDS stops the run in minutes. (The QC node's missing values were
  already counted per chromosome by `GLOWpipeline::prepare_scan_run()`.)
- **`01-run_info.tsv`.** Stage 1 writes the run as one key/value table beside the snapshot: the
  type, the chromosomes, the cohort size and trait, the filter, the inflation factor and where
  it came from (`none`, `manual`, `estimated`, `report`), the type's settings, the models and
  the package versions. The workflow returns it as `run_info`.
- **The example's self-test.** `cloud/terra/run-example-wdl.sh` runs seven workflow runs (the
  gene run calibrated from the report is new), and `check-wdl-example.R` also checks the named
  outputs on disk, the monitor's exit lines, `run_info` and the report, and takes `--tolerance`
  for a reference computed in another environment (the p-value tables stay exact).
- **Documentation.** `cloud/README.md` gains the data-table columns, the outputs per run type,
  the sizing section (workers, memory, the monitor's lines), the defaults to check before a real
  run (`rare_maf_cutoff = 0.5`, the 100 kb windows every 10 kb, `staar_enabled = false`,
  `min_variants = 1`, which follow GLOW's design rather than a rare-variant pipeline's, with the
  lines that switch either route to STAARpipeline's settings) and what Stage 1 checks. The
  workflow's default image is `ghcr.io/zwulab/glowanalyses:v0.3.1`.

# GLOWanalyses 0.3.0

Requires GLOWr 0.3.3 and GLOWpipeline 0.2.1. GLOWr 0.3.3 and GFisher 0.3.2 change the
p-values of the optimal-weight tests (the normal approximation of the weights on the
standardized scale; the HYB method applying unequal weights once), so a scan run with this
release does not reproduce the optimal-weight columns of a scan run with the earlier
packages; the equal-weight tests, the single-variant scan and every statistic are unchanged.
See the NEWS of the two packages.

- **GLOW in the cloud: the official Docker image and the Terra workflow (`cloud/`).**
  `cloud/terra/glow_scan.wdl` runs the scan stages on Terra, or on any machine with miniwdl and
  Docker or Apptainer, as thin tasks over the shipped templates. `region_type = "single_variant"`
  gives the per-variant scan with its calibration (`02-single-variant`), and `gene`, `window`,
  `coding` or `noncoding` give the SNV-set scan with its aggregation (`03-snv-set`, `04-summary`).
  Its inputs carry the configuration fields' names, and its tables equal the SLURM route's on
  the bundled example (`cloud/terra/check-wdl-example.R`, run by `cloud/terra/run-example-wdl.sh`).
  The image `ghcr.io/zwulab/glowanalyses:v0.3.0` is built by `.github/workflows/docker-image.yaml`
  at the release tag from the recipes in `cloud/docker/`. It carries R 4.3.3, Bioconductor 3.18,
  the STAAR family, GFisher 0.3.2, GLOWr 0.3.3, GLOWpipeline 0.2.1 and this tree at
  `/opt/GLOWanalyses`. `.dockstore.yml` registers the workflow on Dockstore. README: `cloud/README.md`
  and the section "Running on Terra". Derived from the Terra workflows written by Alberto
  Brusati (Landers lab, UMass Chan Medical School), 2026.
- **`02-single-variant/compute-ld-scores.R`: per-chromosome units and `--combine`.** `--chr N`
  now writes `ld_chr<N>.rds` beside `ld_scores_path` instead of overwriting the shared table.
  The new `--combine` assembles the per-chromosome files into `ld_scores_path`, the shape
  `marginal-scan.R` already had, so that one machine per chromosome works. The all-chromosome
  call is unchanged, row names included.
- **`ld_ref_maf_min` (`base-config/single_variant_base.R`, default 0).** Above 0, only the SNVs
  with MAF at least this value enter the LD scores, as targets and as neighbors, read from the
  run's single-variant table. On whole-genome data the all-SNV table is out of reach, and
  restricting the reference SNVs can change the estimated LDSC intercept, so the threshold
  travels with the per-chromosome LD files and the calibration should be validated on the
  cohort. The default keeps every SNV, so no existing run changes.
- **The example's PI ensemble also as `data-example/cohort/pi_models.tar.gz`**, the form the
  workflow takes (`build-intermediates.sh` writes it).
- **`run-example.sh` also runs the bundled window and coding example runs**
  (`runs/example/snv-set-window`, `runs/example/snv-set-coding`) after the gene run, so the
  example exercises three of the four set-based region types; `run-example-wdl.sh` runs them
  through the workflow too, plus a coding run with the STAAR comparator, and
  `check-wdl-example.R` compares the tables.
- **Fixes found by the review round before this release.** `03-snv-set/prepare.R` records
  `top_k` in the Stage-1 snapshot (`04-summary/aggregate.R` reads the snapshot, so a `top_k`
  set in `config.R` was ignored in favor of 100 unless `--top-k` was passed); it refuses
  `output_format = "fst"` outside window runs (the gene, coding and noncoding tables are
  written under `.csv` names and read as CSV, so `fst` there failed at aggregation after the
  scan); and with `staar_enabled` off it records `staar_modes` as empty and
  `staar_bakein_glow` as off, so a GLOW-only coding run no longer stops in `run-coding.R` or
  in the aggregator (the manual `staar_bakein_glow <- FALSE` the README asked for is no
  longer needed). `compute-ld-scores.R --combine` refuses a per-chromosome file computed with
  other LD settings than the configuration's, and a MAF threshold that no SNV of a chromosome
  passes stops the script with a message naming the chromosome.

# GLOWanalyses 0.2.0

Requires GLOWr 0.3.0 and GLOWpipeline 0.2.0.

- **The gene-centric non-coding masks, `region_type <- "noncoding"`.** One SNV set per gene
  and mask, the records that the FAVOR annotation links to the gene wherever they lie on the
  chromosome, as in STAARpipeline's gene-centric non-coding analysis: the seven regulatory
  masks of protein-coding genes (`upstream`, `downstream`, `UTR`, `promoter_CAGE`,
  `promoter_DHS`, `enhancer_CAGE`, `enhancer_DHS`) and the `ncRNA` mask, listed in the new
  field `noncoding_categories` of `base-config/snv_set_base.R`. The entry script
  `03-snv-set/run-noncoding.R` runs one chromosome, with `slurm/run-noncoding.sh` as its
  array template, recognized by `slurm/submit-throttled.sh`. It writes
  `results-glow/glow_noncoding_chr<N>.csv` with the mask as `category`. When `staar_enabled`
  is on it carries the STAAR comparator, as the gene run does. `04-summary/aggregate.R`
  writes `glow_noncoding_all.csv`, `noncoding_merged.csv` (a gene's masks together), the two
  significance tables, `noncoding_thresholds.csv` (the Bonferroni cutoff per mask and across
  the configured masks, each with its denominator, STAARpipeline's published constants beside
  them) and `noncoding_completeness.csv` (one row per chromosome with its status and the cells
  attempted, tested, skipped and errored). It aggregates CSV tables only. `plots.R` writes the per-mask QQ plots and lambda
  table. The masks read five aGDS fields that `00-data-prep/annotate-favor.R` writes by
  default: `genecode_comprehensive_category`, `genecode_comprehensive_info`, `genehancer`,
  `cage_tc` and `rdhs`. The two promoter masks need the `TxDb.Hsapiens.UCSC.hg38.knownGene`
  and `GenomicFeatures` packages and an hg38 aGDS. README: "How GLOW forms SNV sets".
- **Size-aware scheduling for runs with `n_cores` above 1.** No field changes. A chromosome's
  regions are now scheduled over the forked workers by their record counts (heavy regions alone,
  the rest in cost-bounded blocks, heaviest first) instead of in blocks of 25 consecutive regions,
  so a task that owns a chromosome and several cores keeps them busy. The written tables are
  unchanged. The per-chromosome log gains one scheduler line per regions table (one per mask
  for the non-coding type) and one line per finished block. On chromosome 22 of the 1000
  Genomes Project with twenty workers on a node allocated to the job alone, the new order took
  55 and 56 percent of the old time in two repeats; with four workers the two were within about
  1 percent of each other.
  A gene list with a repeated label now stops the run at any `n_cores`. README: "Cores and large
  genes".

# GLOWanalyses 0.1.6

Requires GLOWr 0.2.3 and GLOWpipeline 0.1.2, which carry the engine changes behind the two new
groups of fields: the variant index and the faster LD step, the SPA fallback for collapsed
columns, large-gene segmentation and forked workers (see their NEWS files).

- **Large-gene segmentation for gene runs.** Two fields in `base-config/snv_set_base.R`,
  `gene_max_records` (default 10,000) and `gene_segment_bp` (default 10,000), written into the
  Stage-1 snapshot by `03-snv-set/prepare.R`. A gene with more GDS records in its span than the
  first is tiled into segments of the second, each reported as its own row labeled
  `<gene>_seg<k>`. On whole-genome-sequencing data this bounds the per-gene cost and memory;
  genotyping-chip data never reaches the default. Set `gene_max_records <- Inf` to never
  segment.
- **`n_cores`** in `base-config/snv_set_base.R` and `base-config/single_variant_base.R`
  (default 1), passed through the snapshot to `run_scan_unit()` and by
  `02-single-variant/marginal-scan.R` to `marginal_scan()`: the forked workers one task uses
  when it owns a chromosome and several cores, as on a cloud platform. The SLURM templates stay
  at one core per task; raise `--cpus-per-task` and `--mem` together with the field. README:
  "Cores and large genes".

# GLOWanalyses 0.1.5

Requires GLOWr 0.2.0, as 0.1.4 does.

- **Two optional path fields for the FAVOR annotation step.** `favor_input_gds_dir` makes
  `00-data-prep/annotate-favor.R` read a GDS tree other than `<base_name>_gds/`, and
  `favor_output_dir` makes it write one flat tree `<dir>/{gds,csv,provenance,logs}` with no
  `<match>` level, which `compute-pcs.R` and `assemble-pheno-covar.R` then read through the
  shared path resolver instead of rebuilding the path themselves. Left unset, every existing
  config runs exactly as before. The provenance snapshot now records the resolved input and
  output directories. Use the fields to place a re-annotation beside, not over, an earlier
  annotated tree, without a link to the input inside the output folder.
- **The STAAR comparator for gene runs.** `staar_enabled <- TRUE` now works for
  `region_type = "gene"` as it does for window and coding runs: `03-snv-set/prepare.R` fits
  the STAAR null model on the same phenotype and covariates as GLOW's null model, and
  `03-snv-set/run-gene.R` runs `STAAR::STAAR()` on GLOW's own variant set of every gene. The
  per-chromosome tables `results/glow_chr<N>.csv` and the aggregated tables then carry
  `STAAR_O`, `ACAT_O` and the SKAT, Burden and ACAT-V family omnibi at both beta weights
  (`STAAR_S_1_25`, `STAAR_S_1_1`, `STAAR_B_1_25`, `STAAR_B_1_1`, `STAAR_A_1_25`,
  `STAAR_A_1_1`) beside the GLOW tests, weighted by the same annotation features as GLOW's
  variant-importance score when `staar_anno_features` is left `NULL`. `write_staar_detail
  <- TRUE` adds the per-gene STAAR grid as a sidecar under `staar_detail/`. Gene runs of
  earlier versions left these columns empty. STAAR is needed only when the switch is on.
- **README, from user feedback.** The three configuration layers and which ones to edit
  (`base-config/` is overwritten by `git pull`, so your changes go in `_cohort.R` and the
  stage's `config.R`), the two output locations (`data_root` for the `00` data assets,
  `runs/<name>/<stage>/outputs/` for the analysis stages), and the packages to install
  beyond the three GLOW packages, `SNPRelate` for the PCA step among them.

# GLOWanalyses 0.1.4

Requires GLOWr 0.2.0: the annotation templates call arguments of `annotate_favor()` that 0.1.x
does not have.

- **FAVOR 2.0 and the corrected matcher in the annotation templates.**
  `00-data-prep/annotate-favor.R` and `01-training/pi/annotate-pi-case.R` read the FAVOR 2.0
  per-chromosome Parquet files as well as the FAVOR v1 CSV chunks (config field
  `favor_db_format`, default `"auto"`). They record a release label in the aGDS
  (`favor_release`) and read `chroms` as character, so that X and Y pass unchanged. They pass
  the new `favor_rsid_policy` field (`"require"`, the default, or `"record"`) to GLOWr 0.2.0's
  `annotate_favor()`, whose corrected matching records the outcome of every variant
  (`match_tier`) and the rsID evidence (`rsid_check`). The base configs document the fields.
  The README's FAVOR paragraph describes the tiers, the check and the policy.
- **The bundled example** (`data-example/`) is regenerated with GLOWr 0.2.0, so its aGDS
  carries `annotation/info/favor_match_tier`, `annotation/info/favor_rsid_check` and the
  provenance attributes that the corrected annotator writes.

# GLOWanalyses 0.1.3

- **Native STAARpipeline runner templates**, completing the comparison interface
  that `prepare.R` half-shipped (it already fits the native null models for
  coding runs):
  - `03-snv-set/run-staar-coding.R` + `slurm/run-staar-coding.sh`: the native
    `Gene_Centric_Coding` engine (promoted from the validated research driver),
    run once per SPA mode, writing the `*_native` columns the comparison
    aggregation consumes.
  - `03-snv-set/run-staar-window.R` + `slurm/run-staar-window.sh` (new): native
    STAARpipeline sliding-window analysis **unit-matched to the GLOW window
    grid** (the same generator + chunk-ownership rule as the GLOW scan, so rows
    join the GLOW table 1:1 by `label`; STAARpipeline's own window generator is
    not used because it hard-codes a half-window step). Enabled per run by the
    new `staar_native <- TRUE` config field: `prepare.R`'s window branch then
    also fits both native null models and writes the annotation catalog. The
    new `staar_native_use_annotation` field (default TRUE) turns the native
    engine's annotation weighting off for aGDS files with incomplete
    annotation coverage, where the native path's raw (un-imputed) annotation
    reads would otherwise crash the SKAT eigendecomposition.
  STAARpipeline is a documented optional dependency of these templates.

# GLOWanalyses 0.1.2

- The `00-data-prep` steps no longer print `sh: 1: git: not found` at the end of
  the SLURM error log on compute nodes where git is not on `PATH`. The
  provenance git-SHA probe in `00-data-prep/_dataprep_lib.R` now silences the
  probe's stderr and records `Git SHA: unknown` in the output README (was `NA`)
  when git is unavailable. Cosmetic only: analysis outputs were never affected.

# GLOWanalyses 0.1.1

- **SLURM array templates for the two per-chromosome `00-data-prep` steps**:
  `slurm/run-plink-to-gds.sh` and `slurm/run-annotate-favor.sh` (FAVOR annotation
  is the heaviest data-prep step), one task per chromosome over the scripts'
  existing `--chr` entry points. An array task writes its own provenance
  snapshot (`provenance/config_snapshot_chr<N>.rds`) instead of all tasks
  overwriting one.
- **Configurable, robust R-environment activation inside every `slurm/run-*.sh`
  job** (shared `slurm/_job_lib.sh`): `GLOW_CONDA_ENV` (default `r_env`) or
  `GLOW_RSCRIPT` replace the hardcoded `conda activate r_env`. Fixes an abort
  before R started — `...: unbound variable` from conda's compiler-package
  activation hooks under `set -u` — on conda envs that ship compilers.
- Clearer failures: a missing run config reports the path and the working
  directory (was `file.exists(config_path) is not TRUE`); a template run outside
  a SLURM array says so instead of `SLURM_ARRAY_TASK_ID: unbound variable`.
- `00-data-prep/annotate-favor.R`: the documented default `favor_features <- NULL`
  no longer fails with `could not find function .default_favor_features`.
- Documentation: the `00` data-prep-on-SLURM snippet and the R-environment knobs
  in the README; PI wording unified to "variant-importance score".

# GLOWanalyses 0.1.0

Initial public release.

GLOWanalyses is a runnable, cohort-agnostic collection of **workflow templates**
for the GLOW methodology (a scripts directory, not an R package). It is the thin
glue that drives the `GLOWr` (methods) and `GLOWpipeline` (orchestration +
summary) packages across a whole-genome analysis.

- **Five numbered stages**, each a thin `--config` wrapper over the packaged
  methods: `00-data-prep/` (PLINK->GDS, FAVOR annotation to aGDS, PCs,
  pheno/covariate assembly),
    `01-training/` (B effect-size distribution + PI
  variant-importance ensemble),
    `02-single-variant/` (per-variant scan, in-sample LD
  scores, genomic-control calibration),
    `03-snv-set/` (a shared scan for
  the `gene` / `window` / `coding` styles), and
    `04-summary/`
  (region-type-aware aggregation + Manhattan/QQ plots).
- **Run organization** by `runs/<name>/config.R`: a config sources a documented
  `base-config/*_base.R` and names only the cohort paths + what the run varies;
  the stages derive each run's `outputs/` from the config location.
- **A self-contained synthetic example** (`data-example/` + `runs/example/`) that
  runs the whole chain on a tiny tracked cohort with no external data:
  `bash run-example.sh` runs `02`->`03`->`04` over the
  pre-built assets and a schema self-test. The example validates template
  mechanics (every stage emits a schema-valid, non-empty artifact); its outputs
  are statistically meaningless.
- **SLURM array-job templates** (`slurm/`) for the per-chromosome / per-chunk
  Stage-2 scans, with a throttled submitter for capped queues.
- **Optional shared-output storage**, opt-in per run; the templates default to
  portable local outputs.
