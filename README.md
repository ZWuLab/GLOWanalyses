# GLOWanalyses

GLOW (inteGrative anaLysis using Optimized Weights) tests whether a set of genetic variants,
such as the variants of one gene, is associated with a trait. It weights each variant by two
quantities that are learned before the test: B, the effect size expected at the variant's
allele frequency, and PI, a variant-importance score computed from functional annotation. A
GLOW study therefore has several steps. The genotypes are converted and annotated. The B model
and the PI models are trained. Every variant is tested on its own. Then sets of variants are
tested, by gene, by sliding window or by annotation mask. Each step has its own settings and
reads the files of the steps before it. Run by hand, such a study is a series of scripts and
paths that have to agree with each other, on a workstation, on a cluster and in the cloud.

GLOWanalyses is the set of workflow templates that runs these steps on the methods of the R
packages `GLOWr` and `GLOWpipeline`. A study is described by **one configuration file**, in
YAML, that holds only what differs from the shipped defaults. **One command**,
`glow run <step>`, runs a step, in the user's session or as jobs of the SLURM scheduler. The
folder of every output follows from three names in the configuration, so no step is told where
another step wrote. Each unit of work (a part of a step, or one chromosome of it) keeps a
record of the settings and inputs it ran with. A second run of the same command continues
where the first one stopped, and a run that would replace results made with other settings is
refused. GLOWanalyses is a directory of scripts, not an R package, and it does not change the
statistical methods.

The shipped example is a small synthetic cohort on two chromosomes. Its result tables are the
check of the templates: in version 0.4.0 they equal those of the version before, 0.3.6, which
used one configuration script per step, in all 40 result files, and in all 62 after the data
preparation and the training were run again. The example checks the mechanics of the workflow.
Its numbers carry no statistical meaning.

```bash
# first: the R packages of the section "Setup"
git clone https://github.com/ZWuLab/GLOWanalyses.git
cd GLOWanalyses
bash run-example.sh                              # the example: the single-variant and SNV-set steps, then a check

# a study of the user's own: study.yml is a file the user writes (see "The configuration")
./glow run data_prep      --config study.yml     # GDS files, FAVOR annotation, PCs, phenotype bundle
./glow run training       --config study.yml     # the B model and the PI models
./glow run single_variant --config study.yml     # the per-variant scan and its calibration
./glow run snv_set        --config study.yml     # the gene, window or mask scans, with tables and plots
```

On a workstation the command runs the units of a step one after another. On a cluster, the
line `run: {mode: slurm}` in the configuration makes the same command submit one job per
chromosome or chunk. The same steps run on the Terra cloud platform through the workflow in
`cloud/`. A user of version 0.3.6 starts with `MIGRATION.md`.

## Setup

This section says what to install and how to keep a copy up to date.

**The packages.** The templates call three GLOW packages, which are separate repositories and
are installed first, in this order: `GFisher`, `GLOWr` (0.3.6 or later) and `GLOWpipeline`
(0.3.0 or later):

```r
remotes::install_github("ZWuLab/GFisher")
remotes::install_github("ZWuLab/GLOWr")
remotes::install_github("ZWuLab/GLOWpipeline")
```

The templates also attach these Bioconductor and CRAN packages:

| Package | Needed by |
| --- | --- |
| `SeqArray`, `gdsfmt` | every step that reads a GDS file |
| `SNPRelate` | the principal components of the data preparation |
| `arrow` | the FAVOR annotation with a FAVOR 2.0 Parquet database |
| `STAAR` | the STAAR comparator inside the GLOW scan (`snv_set.staar_enabled: true`) |
| `STAARpipeline`, `SeqVarTools` | the native STAARpipeline comparison only (the parts `staar_coding` and `staar_window`) |
| `TxDb.Hsapiens.UCSC.hg38.knownGene`, `GenomicFeatures` | the two promoter masks of a `noncoding` run |

**The copy.** A clone updates in place and states its version in `.release-source`.

```bash
git clone https://github.com/ZWuLab/GLOWanalyses.git
cd GLOWanalyses
cat .release-source        # the version and the source commit of this copy
git pull                   # a later update, or: git checkout v0.4.0 to pin a version
```

A study's own configuration files can sit anywhere, in or outside the clone, and its results
go to the `output_root` that the configuration names. The shipped files are not edited: the
defaults are changed in the user's own configuration, and `git pull` then never meets a
conflict. `git clean` deletes untracked files, which would include a configuration kept in the
clone, so it is not used here.

**Where the commands run.** `./glow` can be called from any directory, and a relative
`--config` is read from the directory of the call. The script of one part, called directly
(`Rscript 03-snv-set/prepare.R --config study.yml`), is run from the GLOWanalyses directory.

## The configuration

This section says what a configuration file holds and how it is checked.

A configuration is a YAML file. It builds on `base-config/base_config.yml`, which lists every
setting with its default and a description, and it holds only the settings that differ. The
shipped example's configuration, `example.yml`, is such a file. A study of the user's own
needs the name of its dataset, the places of its data, and what it changes:

```yaml
# study.yml: every path and name below is a placeholder for the user's own.
dataset: mystudy                        # the first folder under output_root
output_root: /scratch/mylab/glow        # where every result is written (default: runs/ in the clone)

data:
  plink_dir: /data/mystudy/plink        # per-chromosome PLINK files
  plink_pattern: "mystudy_chr{chr}"     # their names without extension; {chr} is the chromosome
  pheno_csv: /data/mystudy/pheno.csv    # one row per sample
  pheno_id_col: IID                     # the column with the sample identifier of the genotypes
  favor_db: /data/favor                 # the FAVOR database on disk

data_prep:
  outcome: {col: status, map: {"control": 0, "case": 1}}   # the outcome column and its recoding
  covariates: {sex: sex, age: age_at_onset}                # name in the model: column of pheno_csv

training:
  b_raw_path: /data/known-snps.csv      # known trait-associated variants, for the B model
  pi_raw_path: /data/known-snps.csv     # and as the cases of the PI models

snv_set:
  region_type: [gene, window]           # the SNV-set scans to run
```

Two inputs are not shipped and come from outside. The FAVOR annotation database is obtained
from the FAVOR project. The table of known trait-associated variants is the user's own. It has
one row per variant, with its position and alleles. From the study that reported the variant
it carries the allele frequency, the effect size or the p-value, and the sample size.
`data-example/cohort/known-snps.csv` shows the columns, and `training.b_column_mapping` and
`training.pi_column_mapping` name them when a table calls them differently.

**The defaults of the raw inputs are the example's.** The settings of the section `data`, and
the two tables of the section `training`, default to the files of the shipped example, so that
the example runs without a configuration of its own. A study under another `dataset` name
that leaves one of them out does not run on the example's data: the run stops and names the
setting. `data_prep.outcome` and `data_prep.covariates` default to the columns of the example's
phenotype table and are set for any other table.

The file has a few names at the top and seven sections.

| Part | What it holds |
| --- | --- |
| `dataset`, `analysis`, `output_root` | the three names that place every output (next section). An analysis is one set of choices made on a dataset, and its default name is `baseline` |
| `data` | the raw inputs: the PLINK files, the phenotype table, the FAVOR database |
| `shared` | what several steps read: the chromosomes, the PI features, the workers per unit, the FAVOR matching |
| `data_prep` | the file-name patterns, the principal components, the outcome and the covariates |
| `training` | the B model and the PI models |
| `single_variant` | the per-variant scan, the LD scores and the calibration |
| `snv_set` | the region types, the variant filter, the tests and the summary |
| `run` | how a step is run: `local` or `slurm`, and the settings of the SLURM mode |

**Checks.** A configuration is read and checked before anything runs, and every problem is
listed with the file it is in. A misspelled setting is refused and the nearest name is offered
(`snv_set.ld_treshold` is not a registered setting; did you mean `ld_threshold`?). A value of
the wrong type is refused with the type that is expected, and a value outside a fixed list is
refused with the list.

**Relative paths.** The steps run in the GLOWanalyses directory, so a relative path of a data
file or folder in a configuration is read from there, wherever the configuration file lies.
The paths of a study are therefore best written in full. Two paths are read differently: a
relative `--config` on the command line is read from the directory of the call, and the file
that an `extends` line names (see "Several analyses of one dataset") is looked for beside the
file that names it.

**Every setting.** `base-config/base_config.yml` is the place to look a setting up: each one
stands there with its default and a description. `base-config/settings-reference.md` gives the
same settings as tables, with the parts of the workflow that read each one and the Terra input
that sets it. Both files are generated from one registry, `base-config/settings-registry.yml`,
which is also what the checks above read.

**What a run used.** Each step folder holds `config-resolved.yml`, written when a unit of the
step runs: every setting the step used with its value, every derived path in full, and which
files of an earlier analysis were reused. It is the record of what the results were made with.

## The steps and their outputs

This section says what each step does, what it reads and where it writes.

`glow run` has four steps. Each step consists of parts, and a part that is split runs once per
chromosome or once per chunk of the genome. One such run of a part is a unit of work: the
command starts units, and the records and the jobs of the next section are kept per unit.

| Step | Parts, in order | Reads | Writes |
| --- | --- | --- | --- |
| `data_prep` | `plink_to_gds`, `annotate_favor` (each per chromosome), `compute_pcs`, `assemble_pheno` | the PLINK files, the FAVOR database, the phenotype table | GDS files, annotated GDS files (aGDS), principal components, the phenotype bundle |
| `training` | `b_prepare`, `b_estimate`, `pi_prepare`, `pi_annotate`, `pi_train`, `pi_evaluate` | the tables of known variants, the FAVOR database, the annotated GDS files | `b_model.rds`, the ensemble `pi_models/`, and an AUC table with a ROC plot |
| `single_variant` | `scan` (per chromosome), `scan_combine`, `ld_scores` (per chromosome), `ld_scores_combine`, `calibrate` | the phenotype bundle, the annotated GDS files | `single_variant_all.csv`, the LD scores, the calibrated table and the calibration report |
| `snv_set` | `prepare`, `scan` (per chromosome, or per chunk for windows), `aggregate`, `plots`, and the two optional native STAARpipeline parts | the phenotype bundle, the annotated GDS files, the B model, the PI models | one table per unit, the aggregated tables, the top hits, Manhattan and QQ plots |

The `snv_set` step runs once for each region type that `snv_set.region_type` lists.

**Where the outputs are.** The folder of every output follows from `output_root`, `dataset` and
`analysis`. The data of a dataset are shared by all its analyses, and everything else belongs
to one analysis:

```
<output_root>/<dataset>/
  data/gds/                          the converted GDS files
  data/gds_favor/<match>/gds/        the annotated GDS files (<match> is shared.favor_match_method)
  data/gds_favor/<match>/csv/        the annotation as tables
  data/pcs/                          the principal components
  <analysis>/
    data_prep/                       the phenotype bundle pheno_covar.rds
    training/b/                      b_model.rds
    training/pi/                     pi_models/ and the evaluation
    single_variant/                  results/, plots/, ld_scores/
    snv_set/<region_type>/           shared/, the unit tables, aggregated/
```

Every step folder also holds `config-resolved.yml`, a folder `records/` with one small file
per finished unit, and its logs. The default `analysis` is `baseline`, and the default
`output_root` is `runs/` in the clone.

**The tables to read first.** The paths below are inside the folder of the analysis,
`<output_root>/<dataset>/<analysis>/`.

| Question | File |
| --- | --- |
| Which variants are associated on their own? | `single_variant/results/single_variant_all_calibrated.csv` |
| How large is the inflation, and how was it corrected? | `single_variant/results/single_variant_calibration.txt`, a report to read. The same facts are in `single_variant_calibration.tsv`, the file that `snv_set.calibration` takes |
| Which genes are associated? | `snv_set/gene/aggregated/glow_results_all.csv` and `glow_top_hits.csv` |
| Which windows, and which loci do they form? | `snv_set/window/aggregated/scan_sig_windows.csv` and `scan_loci.csv` |
| Which coding masks of which genes? | `snv_set/coding/aggregated/coding_merged.csv` |
| Which non-coding masks, at which cutoff? | `snv_set/noncoding/aggregated/noncoding_merged.csv` and `noncoding_thresholds.csv` |
| What was a step run with? | `config-resolved.yml` in the step folder |

## Running a step

This section says what the command does, what a second run does, and how a cluster is used.

```
glow run <step> [--config <file>] [--part <name>] [--force] [--dry-run] [--wait]
```

The command line names the step and the configuration and carries no setting. `--part` runs
one part of the step. `--dry-run` shows every unit of the step with its state and starts none.
`--force` runs again what is finished. Without `--config` the command runs the shipped example.

**A second run.** Every unit writes a completion record when it has finished: the settings it
read, the files it read, the versions of the packages, and the files it wrote. A second run of
the same command compares each unit with its record and skips the units that are complete. So
a step that stopped is continued by running the same command again, and a unit whose result
file was deleted runs again while the others stay.

**Results made with other settings.** A unit whose record states other settings than the
current configuration is not replaced. The command stops, prints each difference with the
recorded and the current value, and names the two ways forward: a new `analysis` name, which
keeps the finished results, or `--force`, which replaces them. A unit that read a part which
was produced again since then is out of date, and it runs again without a question.

**What a step waits for.** A step that needs a part which is not finished stops before any
work and names the step and the part to run first.

**On a SLURM cluster.** With `run.mode: slurm`, the command submits one job per unit that has
to run and returns. A unit that reads another unit of the same step is submitted to start
after that job. A second call leaves alone what is finished or queued and submits what is left.
So the same command is run again until it reports that every unit is finished. `--wait` makes
one call stay, submit the units that become known later (the chunks of a window scan are known
after its preparation), and stop with an error when a job failed. With `--force`, a unit counts
as finished only when its own job wrote its record: the record names the job, so the results of
an earlier run, or of another run in the meantime, do not pass for the work of a job that failed.

```yaml
run:
  mode: slurm
  slurm:
    options: "--partition=short --account=mylab"   # placeholders: the options of the user's own site, added to every sbatch call
    conda_env: r_env                               # the conda environment that a job activates (the default)
    max_queued: 200                                # submit nothing while the user has this many tasks in the queue
    resources:
      default: {time: "12:00:00"}                  # for every part
      annotate_favor: {mem: 64G}                   # for one part
```

A job sets up R by activating the conda environment `run.slurm.conda_env`, which is `r_env`
unless the configuration names another, or by calling the program that `run.slurm.rscript`
names. The job scripts are written to `jobs/` in the step folder and the
logs to `logs/slurm/`, and `--dry-run` writes the scripts without submitting them. A job checks
that the configuration is the one it was submitted with, so a configuration that is edited
while jobs are queued stops those jobs, and the command is run again. Memory, time and cores
are set in `run.slurm.resources`, not in `run.slurm.options`. The command gives the name, the
log, the cores, the memory and the time of a job to `sbatch` itself, and it starts `sbatch`
without the `SBATCH_*` variables of the environment that would set them. Variables of the site,
such as `SBATCH_PARTITION`, are passed on. `--cpus-per-gpu` is refused in `run.slurm.options`,
because `sbatch` does not take it beside the cores that the command asks for. `--time-min` is
accepted: with it the scheduler may start a job earlier under a shorter time limit, never a
longer one than `run.slurm.resources` gives. `max_queued` is read when a call looks at the queue:
a call submits nothing while the user has that many tasks queued. It is not a hard limit when
two calls for two datasets start in the same moment.

## Several analyses of one dataset

This section says how a second analysis reuses the work of the first.

An analysis is one set of choices on a dataset. A second analysis is a second configuration
file that names the first in its `extends` line, gives itself a new `analysis` name and states
only what it changes. The file named by `extends` is looked for beside the file that names it:

```yaml
# strict.yml: builds on study.yml and changes one setting
extends: study.yml
analysis: strict
snv_set:
  ld_threshold: 0.8
```

Its results go to `<output_root>/<dataset>/strict/`. It reuses, without copying, every part of
the analysis it extends that does not depend on what it changed. Here those are the phenotype
bundle, the two sets of models and the single-variant results, while the SNV-set step runs
again. When the change reaches an earlier step (other PI features need other PI models), the
command stops at the step that would read the old result and names the step to run first. A
step folder is reused whole or not at all, and only from an analysis that finished it.

## Inputs that already exist

This section says how a study starts from files it already has.

A study does not have to run every step. A setting that names a finished file takes the place
of the part that would produce it, and that part is left out of the run.

| The study already has | Setting |
| --- | --- |
| GDS files, not yet annotated | `data_prep.favor_input_gds_dir`, with the pattern `data_prep.gds_pattern` |
| annotated GDS files | `shared.analysis_gds_dir` and `shared.analysis_gds_pattern` |
| a table of principal components | `data_prep.pcs.path` |
| a phenotype bundle | `shared.pheno_path` |
| a B model, or a B function | `snv_set.b_model_path`, or `snv_set.b_func` |
| trained PI models | `snv_set.pi_model_dir` |
| LD scores | `single_variant.ld_scores_path` |

A study with annotated GDS files, a phenotype bundle, a B model and PI models runs
`glow run single_variant` and `glow run snv_set` and nothing before them.

A study that supplies annotated GDS files and still trains its PI models also names the folder
from which the control variants are drawn, `training.pi_control_source`. By default that is the
folder the annotation of the dataset writes, which a run with supplied files does not have. The
folder of the supplied GDS files serves when its files are named `chr<N>...gds`, with
`training.pi_control_format: gds`. Without the setting, `glow run training` stops and says so.

## The example

This section says what the shipped example is and what its check covers.

`data-example/` holds a small synthetic cohort on chromosomes 21 and 22, with a small FAVOR
database, a phenotype table and a table of known variants, all tracked. It also holds the files
that the data preparation and the training produce from them (the annotated GDS files, the
phenotype bundle, the B model and five PI models), so the example starts at the analysis steps.

```bash
bash run-example.sh                  # the single-variant step and the gene, window and coding scans, then the check
REBUILD=1 bash run-example.sh        # also regenerate the cohort and run the data preparation and the training
DRYRUN=1  bash run-example.sh        # show every unit and its state, and start none
```

The results are written to `runs/example/baseline/`. The last step,
`data-example/check-example.R`, checks that each stage wrote its tables with the expected
columns and at least one row, and that every unit of the steps that ran is complete by its
record. It does not judge p-values. A second `bash run-example.sh` finds every unit complete
and does no work. `data-example/README.md` describes the cohort and how its
pre-built files are renewed.

## Running on Terra

This section says how the same steps run in the cloud.

The single-variant scan and the SNV-set scans also run on Terra, or on any WDL engine, without
an installation. `cloud/terra/glow_scan.wdl` runs the same templates inside the image
`ghcr.io/zwulab/glowanalyses:<version>`, which carries R, the three GLOW packages and this
tree at matching versions. A Terra user sets workflow inputs, not a configuration file, and the
inputs carry the names of the settings. `cloud/README.md` has the inputs, the sizing, the
defaults to check before a real run, and what changed for a Terra user in this version.

## How GLOW forms SNV sets

This section says which variants enter a set under each region type.

`snv_set.region_type` lists the region types to run, and one shared scan
(`GLOWpipeline::run_scan_unit()`) serves all four. A scan unit is one chromosome (for the
`gene`, `coding` and `noncoding` types) or one genome chunk (for `window`). A unit holds many
regions, and GLOW tests one SNV set per region. The set is formed in two steps that the code
keeps apart. The **region** says which variants are candidates. The **filter**
(`snv_set.filter_spec`: the QC label, the variant type, the MAF and MAC cutoffs, the minimum
set size) says which candidates stay. The four types differ in how the region is defined:

| `region_type` | one set per | candidates | how a variant is tied to the set | unit | table of a unit |
|---|---|---|---|---|---|
| `gene` | gene, or segment of a gene above `gene_max_records` | every record in the gene's or segment's span | by position | chromosome | `results/glow_chr<N>.csv` |
| `window` | sliding window | every record in the window | by position | chunk | `results/scan_chunk_<K>.csv` |
| `coding` | gene and coding mask | the records in the gene's span whose coding annotation matches the mask | by position, then by annotation | chromosome | `results-glow/glow_coding_chr<N>.csv` |
| `noncoding` | gene and non-coding mask | the records anywhere on the chromosome whose annotation names the gene and matches the mask | by annotation only | chromosome | `results-glow/glow_noncoding_chr<N>.csv` |

**`gene`.** One set per gene of the bundled protein-coding gene table (18,445 genes). On
whole-genome sequencing data a gene with more records than `snv_set.gene_max_records` is tiled
into segments of `snv_set.gene_segment_bp`, each tested as its own row labeled
`<gene>_seg<k>` (see "Cores and large genes").

**`window`.** One set per sliding window of `snv_set.window_size` base pairs every
`snv_set.step_size`, over the whole chromosome, so non-coding regions are covered by position
without any mask.

**`coding`.** One set per gene and coding mask, the records in the gene's span whose coding
annotation matches the mask. The built-in masks follow STAARpipeline: `plof`, `plof_ds`,
`missense`, `disruptive_missense`, `synonymous`, `ptv`, `ptv_ds` (`snv_set.categories`).

**`noncoding`: the gene-centric non-coding masks.** One set per gene and mask, the records
that the annotation links to the gene wherever they lie on the chromosome, as in
STAARpipeline's gene-centric non-coding analysis. A variant upstream of a gene lies outside
the gene body, so gene coordinates cannot locate it. The link comes from the FAVOR annotation
of the record. The masks use these terms:

- GENCODE: the gene annotation whose category (upstream, downstream, UTR, exonic, intronic
  and so on) and gene names FAVOR carries per record.
- UTR: an untranslated region of a transcript. ncRNA: a non-coding RNA gene.
- CAGE and DHS: experimental marks of active promoters and of open chromatin, recorded by
  FAVOR as the `cage_tc` and `rdhs` fields.
- GeneHancer: a database of enhancer elements with the genes they are connected to,
  recorded as the `genehancer` field.
- TSS: a transcription start site.

The eight masks (`snv_set.noncoding_categories`):

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
`ncRNA`. Under the promoter and enhancer masks it joins one gene's set.

The five aGDS fields the masks read are `genecode_comprehensive_category`,
`genecode_comprehensive_info`, `genehancer`, `cage_tc` and `rdhs`, which the data preparation
writes by default. Both FAVOR layouts are read, the FAVOR v1 strings and the FAVOR 2.0
serialization of the annotator. The FAVOR 2.0 layout does not keep the group boundary of a
record with two categories, and such a record is then linked to every gene it names under each
of its masks. The two promoter masks take the transcription start sites from
`TxDb.Hsapiens.UCSC.hg38.knownGene` (with `GenomicFeatures`), genome-wide, so a variant can lie
in gene A's window and be linked to gene B by its annotation. The aGDS has to be on hg38 for
these two masks. We call the eight masks the gene-centric non-coding masks, following
STAARpipeline. Four of them (`upstream`, `downstream`, `UTR`, `ncRNA`) select by GENCODE
category. The promoter and enhancer masks select by the CAGE, DHS and GeneHancer annotations
instead, so they can include exonic variants. The run's `filter_spec` applies to every mask
without an added annotation clause, and the sets are LD-pruned, collapsed and tested exactly as
the other types.

GLOW reads the annotation by its meaning rather than by STAARpipeline's string operations,
and a set can differ from STAARpipeline's in four cases. A record annotated
`upstream;downstream` joins both masks. A UTR record naming two genes joins both genes'
sets. A hyphenated gene name is kept whole in the promoter masks. An ncRNA record with
several gene groups links them all. We validated this on chromosome 22 of the 1000 Genomes
Project (3,202 samples, 427 protein-coding and 483 ncRNA genes). The sets formed with
STAARpipeline 0.9.8's own string operations (`rule = "staar"` of
`GLOWr::define_regions_noncoding()`) matched the native selection variant for variant, and
the default rule differed only in those four cases. The table of a unit has the coding
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

## Settings in detail

This section explains the settings whose meaning a default does not tell: the FAVOR
annotation, the SNV-set scan, and the use of several cores.

### The FAVOR annotation

`data.favor_db` is a FAVOR v1 CSV directory (the `chr*_*.csv` chunks) or a FAVOR 2.0 directory
(one `chromosome_<N>.parquet` per chromosome). The FAVOR 2.0 files are obtained from the FAVOR
team. `shared.favor_db_format` is `auto` by default and detects which one it is, and the
Parquet format needs the R package `arrow`. `shared.favor_match_method: exact` looks each
variant up by its own CHR-POS-REF-ALT. `flexible`, the default, also tries, for SNVs, the key
with REF and ALT swapped and the key on the opposite strand, which genotyping-chip data
converted through PLINK can need. The genotypes and alleles are never changed. Only the lookup
key is. Every variant gets an outcome, `match_tier` in the CSV and
`annotation/info/favor_match_tier` in the aGDS: `exact`, `swapped`, `flipped`,
`flipped_swapped`, `position_only` (a key without alleles), or a reason it was not annotated
(`rsid_conflict`, `unmatched_alt`, `unmatched_ref`, `uncovered`, `no_database`,
`unsupported_allele`). A swapped or flipped outcome says which lookup found the row. It is not a
check of the strand. The evidence for a match is the rsID. FAVOR carries dbSNP's rsID on the
specific REF-ALT row dbSNP knows. So every matched variant also records whether the chip's rsID
(the GDS `annotation/id`) equals the FAVOR row's, in `rsid_check` (`same`, `differs`,
`chip_none`, `favor_none`) and `annotation/info/favor_rsid_check`. `shared.favor_rsid_policy`
is `require` by default. Under it, a swapped or flipped match that the rsID does not confirm is
withheld as `rsid_conflict`. So a rewritten key is accepted only when dbSNP confirms it or the
data carries no rsID to check. `record` keeps every match and only records the check, for
diagnostics. Position-only keys are annotated with the average over the SNV rows at the
position. `shared.favor_release` is a label recorded in the aGDS, since the FAVOR files do not
state their release. The vignette `vignette("favor-annotation", package = "GLOWr")` explains
the matching rule, the tiers and the rsID evidence. The annotation is the heaviest part of the
data preparation, about 32 GB at its peak and 15 minutes per chromosome measured against the
FAVOR Essential database, which is why its job asks for 40 GB by default.

### The SNV-set scan

The settings below are in the section `snv_set` unless another section is named.

- **The weights.** The scan needs the PI models (`pi_model_dir`, by default those of the
  training step) and one B source: the B model of the training step (`b_model_path`), or a
  function of the allele frequency written in the configuration (`b_func`).
- **The variant filter.** `filter_spec` holds `rare_maf_cutoff`, `variant_type`, `min_mac` and
  `min_variants`. The shipped values (0.5, `SNV`, 1 and 1) keep the whole frequency spectrum
  and one-variant sets, which is GLOW's design and not STAARpipeline's rare-variant convention.
  `cloud/README.md`, "Defaults to check before a real run", has the lines that set the
  rare-variant convention.
- **Preprocessing.** `ld_threshold` (default 0.9): LD pruning removes one variant of every pair
  whose absolute correlation is above it. `mac_threshold` (default 10, 0 disables): variants
  with a cohort MAC strictly below it are collapsed into one unit. `collapse_method` (`mean`,
  `max` or `sum`, default `mean`): how the B and PI weights of a collapsed unit are aggregated,
  while its genotypes are summed. The three defaults are those of
  `GLOWr::glow_preprocessing_defaults()`, and `shared/01-run_info.tsv` of a run records the
  values it applied.
- **Calibration.** `calibration` is off by default. It takes the calibration report of a
  single-variant run (`single_variant_calibration.tsv`), or the two files that report is
  computed from. The preparation refuses a report of another phenotype bundle.
- **The STAAR comparator.** With `staar_enabled: true` the preparation fits the STAAR null
  model on GLOW's phenotype and covariates, and every row of every region type also carries
  `STAAR::STAAR()` on GLOW's own variant set: `STAAR_O`, `ACAT_O` and the six family omnibus
  tests, named `STAAR_O_glowG` and `ACAT_O_glowG` in the coding aggregates. The comparator
  hands STAAR the genotype matrix as the GDS file stores it, and STAAR flips it to the minor
  allele inside. STAARpipeline flips before its call, so on an identical variant set the two
  STAAR-O values agree to about 1e-8 relative, not to machine precision. `staar_anno_features`
  names its annotation-weight columns (null takes the PI features).
- **The native STAARpipeline comparison.** With `staar_enabled: true`, a coding run also runs
  the part `staar_coding`, STAARpipeline's own engine for coding masks, once per mode of
  `staar_modes` (`spa` and `nospa`; an empty list leaves the part out). A window run does the
  same with the part `staar_window` when `staar_native: true` is set as well. Both parts need
  the packages `STAARpipeline` and `SeqVarTools`. The annotated GDS files need the four fields
  of STAARpipeline's annotation catalog (`genecode_comprehensive_category`,
  `genecode_comprehensive_exonic_category`, `genecode_comprehensive_info` and `metasvm_pred`),
  which the data preparation writes by default. The shipped example annotates a shorter list
  without `genecode_comprehensive_info`, so the two parts do not run on it.
- **Outputs.** `output_format` is `csv`. `fst` is accepted for window runs only. The
  `write_*` settings switch the side files of a scan on. `top_k` is the number of rows of
  `glow_top_hits.csv` in a gene run (default 100). `primary_test`, `primary_test_glow` and
  `alpha` decide what the summaries call significant. `debug_errors: true` puts the failing
  call and the stack of function names into the `message` of an error row.
- **Small test runs.** `max_genes` and `max_windows` cap the regions of each unit, and
  `shared.chroms` limits the chromosomes. `allow_incomplete: true` lets the aggregation run on
  the units that are finished.
- **Window runs.** `window_size`, `step_size`, `chunk_strategy`, `chunk_size_mb`, `merge_gap`
  and `cmac_cutoff`.
- **Non-coding runs.** `noncoding_categories` lists the masks (default all eight). The
  aggregation writes `glow_noncoding_all.csv`, `noncoding_merged.csv`, two significance tables,
  `noncoding_thresholds.csv` and `noncoding_completeness.csv`. The table of a unit has the
  columns `label` (the gene), `chr`, `category` (the mask), `start`, `end` and `midpoint` of
  the candidate span, and the count chain of every GLOW table: `n_variants` (the candidates in
  the file), `n_after_annotation` (after the QC and type settings), `n_after_filter` (after the
  MAF and MAC cutoffs), `n_after_ld` and `n_after_collapse`, `cMAC` and `is_single`. Then come
  the GLOW tests, `GLOW_Omni` the primary one. `noncoding_thresholds.csv` has one row per
  method and mask and one per method for `all`, the configured masks together: `n_tested`, the
  computed `cutoff`, its `cutoff_rule` and the reference constant with its derivation.
  `noncoding_completeness.csv` has one row per requested chromosome: its `status`, the cells
  and genes, the tested and missing primary p-values, and the cells attempted, ok, skipped and
  errored.

### Cores and large genes

This part says how a scan uses several cores, and what happens to a gene with very many
variants. By default each unit of a scan is one process on one core, and on a SLURM cluster one job. On
a cloud platform such as Terra, where a task is a virtual machine that first copies its inputs,
the economical layout is one task per chromosome with the machine's cores filled from inside
R. `shared.n_cores` sets that number. `run_scan_unit()` then splits the chromosome's regions
into blocks and runs them in forked workers, each with its own GDS handle. The single-variant
scan does the same with its variant chunks. The written tables are identical to
the single-core ones. The memory of a task is sized as `n_cores` times the per-region peak, and
`OPENBLAS_NUM_THREADS=1` keeps the workers from oversubscribing a threaded BLAS. For scale: the
`noncoding` unit of chromosome 22 of the 1000 Genomes Project (3,202 samples, eight masks, the
STAAR comparator on, four workers) had a peak R heap of 1.3 GB in the parent process and a
job-wide peak of 8.4 GB. On a 16 GB node, two workers are a safe start, and `n_cores` is raised
once a large chromosome has run. In the SLURM mode a job of a scan asks for `shared.n_cores`
cores by itself, and its memory is raised in `run.slurm.resources`. Forking is not available on
Windows.

The regions of a chromosome are not handed to the workers in genomic order. Genes with many
records sit next to each other, and a block of neighbors would keep one worker busy after the
others have finished. Instead, `run_scan_unit()` reads every region's record count from the
file's index and takes the squared count as the region's cost. A region whose cost reaches the
chromosome's total cost divided by 8 times `n_cores` runs in a block of its own. The remaining
regions are packed, in decreasing cost, into blocks of at most 25 regions whose summed cost
stays within that threshold, and the blocks are dispatched heaviest first. The aim is a wall
time close to the total work divided by the cores. For the gene, window and coding types the
record count is read from the file's index over the region's span. For the non-coding masks it
is the region's candidate count. The budget and the blocks are formed per regions table: once
per unit for a gene or window run, once per chromosome for the coding categories, which share
the gene table, and once per mask for the non-coding type. On chromosome 22 of the 1000 Genomes
Project (3,202 samples) with twenty workers on a node allocated to the job alone, this order
took 55 and 56 percent of the time of the genomic-order blocks in two repeats. With four
workers the two orders were within about 1 percent of each other. The log of a unit shows one
scheduler line per regions table (regions, blocks, how many run alone, the heaviest region) and
one line per block with its region count, its heaviest record count and its wall time, written
after all of the table's blocks have returned, in dispatch order, so the log shows where the
time went. A gene list with a repeated label stops the run, at any `n_cores`, because the
output files of a region are keyed by its label.

On whole-genome sequencing data a gene body can hold tens of thousands of variant records, and
the per-gene cost and memory grow with the square of the variant count. A gene run therefore
tiles every gene with more GDS records in its span than `snv_set.gene_max_records`. The
records are counted from the file's positions, before any filter. The segments do not overlap
and are `snv_set.gene_segment_bp` base pairs long, the last one ending at the gene's end. Each
segment is tested and reported as its own row labeled `<gene>_seg<k>`, with no gene-level
combination. This is how the GLOW methodology paper handled
genes above 1,000 variants, with 10 kb windows. The defaults, 10,000 records and 10 kb, leave
genotyping-chip data untouched. A lower `gene_max_records` suits a large cohort (the genotype
matrix of one gene is 8 bytes times the sample count times the variant count).

## Coming from release 0.3.6

Release 0.3.6 was configured by `config.R` files in three layers and run script by script.
`MIGRATION.md` says where each field, option, script and output folder of that release went,
with one worked example. A `config.R` given to the command is refused with a pointer to that
file.


> Developed with AI assistance; see [AI-USE.md](AI-USE.md).
