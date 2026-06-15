# GLOWanalyses 0.1.0

Initial public release.

GLOWanalyses is a runnable, cohort-agnostic collection of **workflow templates**
for the GLOW methodology (a scripts directory, not an R package). It is the thin
glue that drives the `GLOWr` (methods) and `GLOWpipeline` (orchestration +
summary) packages across a whole-genome analysis.

- **Five numbered stages**, each a thin `--config` wrapper over the packaged
  methods: `00-data-prep/` (PLINK->GDS, FAVOR annotation to aGDS, PCs,
  pheno/covariate assembly), `01-training/` (B effect-size distribution + PI
  pathogenicity ensemble), `02-single-variant/` (per-variant scan, in-sample LD
  scores, genomic-control calibration), `03-snv-set/` (a shared scan for
  the `gene` / `window` / `coding` styles), and `04-summary/`
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
