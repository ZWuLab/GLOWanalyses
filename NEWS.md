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
