# Moving a study from GLOWanalyses 0.3.6 to 0.4.0

GLOWanalyses 0.3.6 was configured by R scripts in three layers and run script by script.
Version 0.4.0 reads one YAML configuration and runs a stage with one command. The statistical
methods are the same, and on the shipped example the two versions give the same result tables.
This guide says where everything of a study on 0.3.6 goes. A `config.R` is no longer read: a
script or the command that gets one stops and points here. Version 0.4.0 needs GLOWr 0.3.6 or
later and GLOWpipeline 0.3.0 or later. A user of the Terra workflow has nothing to move, and
section 11 says what differs there.

## 1. What changed

This section gives the change in one table.

| | Release 0.3.6 | Release 0.4.0 |
| --- | --- | --- |
| Configuration | `base-config/*_base.R`, a `_cohort.R` and one `config.R` per stage and region type | one YAML file that builds on `base-config/base_config.yml` |
| A second analysis | a copy of the run folder with edited files | a second YAML file with `extends` and a new `analysis` name |
| Output folders | beside each `config.R` (`outputs/`), and a `data_root` tree for the data | derived from three names: `<output_root>/<dataset>/data/` for the data of a dataset, and `<output_root>/<dataset>/<analysis>/<step>/` for the rest (section 5) |
| Running a stage | one `Rscript` call per script, in the right order | `glow run <step>` |
| A cluster | job templates in `slurm/`, submitted by hand as arrays | `run.mode: slurm`, and the same command submits the jobs |
| A second run | ran everything again | skips what is finished and refuses to replace results made with other settings |
| Paths between stages | written by hand in each `config.R` | found from the layout |

A result of 0.3.6 is not read by 0.4.0 as a finished step, because it has no completion record.
So there are two ways to go on. A study that is under way is best finished on the version it
started with: `git checkout v0.3.6` in the clone keeps that version, with its base files and
job templates. A new analysis on 0.4.0 can start from the large files of the old study (the
annotated genotypes, the phenotype bundle, the models), by the settings of section 8.

## 2. The move, step by step

This section gives the order of the work. The sections after it have the details.

1. **Update.** Install GLOWr 0.3.6 and GLOWpipeline 0.3.0, then `git pull` in the clone. The
   update removes the shipped `base-config/*_base.R` files, the `config.R` files of the example
   and the job templates of `slurm/`. The user's own `runs/<name>/` folders, with their
   `config.R` files and `outputs/`, are not tracked and stay as they are. Version 0.4.0 ignores
   the whole folder `runs/`, so they no longer show in `git status`.
2. **Write one file.** Start a YAML file with the name of the dataset, the place of the
   results and the raw inputs, which `_cohort.R` and the data-preparation `config.R` held:

   ```yaml
   # study.yml: every path and name is a placeholder for the user's own
   dataset: mystudy                      # was base_name
   output_root: /scratch/mylab/glow      # where the results go; was data_root and the runs/ folders
   data:
     plink_dir: /data/mystudy/plink      # was <data_root>/<base_name>_plink
     plink_pattern: "mystudy_chr{chr}"
     pheno_csv: /data/mystudy/pheno.csv
     pheno_id_col: IID
     favor_db: /data/favor
   ```

3. **Move the settings.** For every field that `_cohort.R` or a stage's `config.R` set, look up
   its setting in section 10 and write it under its section, in the YAML form of section 4.
   The paths from one stage to the next are left out (section 3).
4. **Name the files that already exist**, if the new analysis is to start from them
   (section 8).
5. **Look before running.** `./glow run <step> --config study.yml --dry-run` reads and checks
   the file and shows every unit of the step with its state. A misspelled setting or a value
   of the wrong type is reported here, with the name of the setting.
6. **Run**, with the commands of section 6.

On the shipped example the two versions give identical result tables. For a study, the same
kind of check is quick: run one part on both versions, for example the gene scan of one
chromosome, and compare the two tables, which have the same columns.

## 3. From three layers to one file

This section shows how the files of a study become one.

In 0.3.6 a setting could stand in a base file, in `_cohort.R` or in a stage's `config.R`, and
the last one read won. In 0.4.0 every default stands in `base-config/base_config.yml`, which is
not edited, and the study's file holds what differs. The fields of `_cohort.R` and of all stage
files go into that one file, each under the section that section 10 names for it. A field whose
value is the default is left out.

The window run of the shipped example was these two files in 0.3.6, shortened here:

```r
# runs/example/_cohort.R
data_root <- "data-example/cohort"
base_name <- "example"
chroms    <- c(21L, 22L)
favor_gds_pattern <- "chr{chr}_favor.gds"

# runs/example/snv-set-window/config.R
source("base-config/snv_set_base.R", local = TRUE)
source("runs/example/_cohort.R",     local = TRUE)
region_type  <- "window"
gds_dir      <- file.path(data_root, paste0(base_name, "_gds_favor"), "flexible", "gds")
gds_pattern  <- favor_gds_pattern
pheno_path   <- file.path(data_root, "pheno_covar.rds")
b_model_path <- file.path(data_root, "b_model.rds")
pi_model_dir <- file.path(data_root, "pi_models")
window_size  <- 50000L
step_size    <- 25000L
```

In 0.4.0 the same run is part of `example.yml`:

```yaml
shared:
  chroms: [21, 22]
data_prep:
  favor_gds_pattern: "chr{chr}_favor.gds"
snv_set:
  region_type: [gene, window, coding]   # the three runs of the example in one list
  window_size: 50000
  step_size: 25000
```

Three things went away. The two `source()` lines have no counterpart: a file builds on the base
unless its `extends` line names another file. The five paths have none either: the annotated
GDS files, the phenotype bundle and the models are found from the outputs of the steps that
write them. And the three stage files of the three region types are one list,
`snv_set.region_type`. A setting that only one region type reads, such as `window_size`, stands
beside the others and is ignored by the runs that do not read it.

One case needs a choice. A setting that several region types read has one value per file. If
the old stage files gave it different values on purpose, for example another `ld_threshold`
for the window run than for the gene run, each value needs a file of its own: a second file
that names the first in its `extends` line, has its own `analysis` name, lists its region
type and sets the value (the README, "Several analyses of one dataset"). The example of
0.3.6 set `top_k` to 20 for the gene run and left it at 100 elsewhere. Only the gene run reads
`top_k`, so the one value 20 in `example.yml` is the same example.

## 4. An R value in YAML

This section gives the YAML form of the R values that a `config.R` used.

| In R | In YAML | Note |
| --- | --- | --- |
| `TRUE`, `FALSE` | `true`, `false` | |
| `NULL` | `~` | or the word `null` |
| `10L`, `0.05`, `2.5e-6` | `10`, `0.05`, `2.5e-6` | no `L` |
| `"SNV"` | `SNV` | quotes are needed only around text that YAML would read otherwise, such as `"N"` or `"chr{chr}"` |
| `1:22` | `1-22` | a range, in a list of whole numbers such as `shared.chroms` |
| `c(21L, 22L)` | `[21, 22]` | |
| `c("plof", "missense")` | `[plof, missense]` | |
| `list(sex = "sex", age = "age")` | `{sex: sex, age: age}` | or one line per entry, indented |
| `file.path(root, "pheno.csv")` | `/data/study/pheno.csv` | a path is written out; YAML computes nothing |
| `function(maf) dbeta(maf, 1, 25)` | `"function(maf) dbeta(maf, 1, 25)"` | `snv_set.b_func` takes the text of the function, which may use base R, `stats` and `GLOWr` |

A value is read by the type of its setting when the file is read. So `"0.9"` in quotes is the
number 0.9 for a numeric setting, and `"2"` the integer 2. A value that cannot be read as
that type, or a word that is not one of the allowed choices, stops the run with the name of
the setting.

## 5. Where the outputs went

This section maps the folders of 0.3.6 to those of 0.4.0.

| Release 0.3.6 | Release 0.4.0 |
| --- | --- |
| `<data_root>/<base_name>_gds/` | `<output_root>/<dataset>/data/gds/` |
| `<data_root>/<base_name>_gds_favor/<match>/{gds,csv}/` | `<output_root>/<dataset>/data/gds_favor/<match>/{gds,csv}/` |
| `<data_root>/<base_name>_pcs/` | `<output_root>/<dataset>/data/pcs/` |
| `<data_root>/<base_name>_pheno/<run>_pheno_covar.rds` | `<output_root>/<dataset>/<analysis>/data_prep/pheno_covar.rds` |
| `runs/<name>/training/outputs/` | `<output_root>/<dataset>/<analysis>/training/b/` and `training/pi/` |
| `runs/<name>/single-variant/outputs/` | `<output_root>/<dataset>/<analysis>/single_variant/` |
| the shared LD table at `ld_scores_path` | `<analysis>/single_variant/ld_scores/ld_scores.rds`, unless `single_variant.ld_scores_path` names another place |
| `runs/<name>/snv-set*/outputs/` | `<output_root>/<dataset>/<analysis>/snv_set/<region_type>/` |

Inside a step folder the files keep their names, with one exception. The tables and the plot
of the single-variant scan are `single_variant_chr<N>.csv`, `single_variant_all.csv` and
`single_variant_plots.pdf`. They were `marginal_chr<N>.csv`, `marginal_all.csv` and
`marginal_plots.pdf`.

## 6. The commands

This section gives the command that replaces each script and job template.

| Release 0.3.6 | Release 0.4.0 |
| --- | --- |
| `Rscript 00-data-prep/<script>.R --config <run>/data-prep/config.R`, four scripts | `./glow run data_prep --config study.yml` |
| `Rscript 01-training/b/...` and `01-training/pi/...`, six scripts | `./glow run training --config study.yml` |
| `Rscript 02-single-variant/...`, three scripts | `./glow run single_variant --config study.yml` |
| `Rscript 03-snv-set/prepare.R`, the `run-*.R` scripts and `04-summary/*.R`, per region type | `./glow run snv_set --config study.yml` |
| `sbatch --array=1-22 slurm/run-gene.sh <config>` and the other templates | `run: {mode: slurm}` in the configuration, and the same four commands |
| `CAP=90 bash slurm/submit-throttled.sh <config>` | `run.slurm.max_queued: 90` |
| `GLOW_CONDA_ENV`, `GLOW_RSCRIPT` | `run.slurm.conda_env`, `run.slurm.rscript` |
| `sbatch` flags on the command line | `run.slurm.options` for the options of the site, `run.slurm.resources` for memory, time and cores |
| one script of a stage, for all chromosomes | `./glow run <step> --config study.yml --part <part>` |
| one script, for one chromosome | the script itself, with `--config study.yml --chr <N>` |

The stage scripts are still there, and the command calls them. A script called directly takes
`--config <file>.yml`, the selector of its unit (`--chr`, `--chunk`, `--combine`) and, for an
SNV-set script, `--region-type`.

## 7. Options that became settings

This section lists the values that were given on the command line of a script.

A value that changes what a run computes is a setting in 0.4.0, so that it is in the record of
the run. A script that gets one of these options stops and names the setting.

| Option of 0.3.6 | Setting of 0.4.0 |
| --- | --- |
| `--method` of `calibrate.R` | `single_variant.calibration_method` |
| `--top-k` of `aggregate.R` | `snv_set.top_k` |
| `--allow-incomplete` of `aggregate.R` | `snv_set.allow_incomplete` |
| `--drilldown` of `plots.R` | `snv_set.drilldown` |
| `--max-genes` of the gene, coding and non-coding scans | `snv_set.max_genes` |
| `--max-windows` of the window scans | `snv_set.max_windows` |
| `--n-cores` of the scans | `shared.n_cores` |

## 8. Starting from the files of a study on 0.3.6

This section says how a new analysis reads the large files that an old study already has.

The files of a study on 0.3.6 can be named in the configuration, and the step that would
produce each one is then left out.

| File of the old study | Setting |
| --- | --- |
| the annotated GDS files, `<base_name>_gds_favor/<match>/gds/` | `shared.analysis_gds_dir`, with `shared.analysis_gds_pattern` |
| the converted GDS files, `<base_name>_gds/` | `data_prep.favor_input_gds_dir`, with `data_prep.gds_pattern` |
| the table of principal components | `data_prep.pcs.path` |
| the phenotype bundle | `shared.pheno_path` |
| `b_model.rds` | `snv_set.b_model_path` |
| the folder `pi_models/` | `snv_set.pi_model_dir` |
| the LD table | `single_variant.ld_scores_path` |

## 9. Names that changed, and fields that are gone

This section lists the fields whose name is not simply moved under a section.

| Field of 0.3.6 | In 0.4.0 |
| --- | --- |
| `gds_dir`, `gds_pattern` of the analysis stages | `shared.analysis_gds_dir`, `shared.analysis_gds_pattern` |
| `use_SPA` of the single-variant stage | `single_variant.use_spa` |
| `data_root`, `base_name` | gone: `dataset` and `output_root` place the outputs, and the section `data` names the raw inputs |
| `output_dir` | gone: the output folder is derived |
| `favor_output_dir` | gone: the annotated files are written to `data/gds_favor/<match>/` of the dataset |
| `symlinked_shared_root` | gone: `output_root` says where the outputs go, and can point at shared storage |
| `qc_filters$filter_indels` of the PI case table | gone with GLOWr 0.3.6: the option was never applied |

A field that several stage files set, such as `chroms`, `pi_features` or `favor_match_method`,
is one setting in the section `shared`.

## 10. Every field

This section is the table of every field of 0.3.6 with its setting in 0.4.0: first the fields
that a base file set, in the order of the four files, then the fields that a base file only
described and a run had to set, then two fields that no base file held. A field with an empty
note moved under its section and kept its name and its default. Two notes recur. "Was unset"
means that the base file gave no value and every run set one. "The default is the file of the
shipped example" means that the setting now has a default, which is right for the example
only, so a study sets it.

| Base file of 0.3.6 | Field | Setting of 0.4.0 | Note |
| --- | --- | --- | --- |
| `data_prep_base.R` | `symlinked_shared_root` | none | gone: `output_root` says where the outputs go |
| `data_prep_base.R` | `chroms` | `shared.chroms` | one setting for the steps that share it |
| `data_prep_base.R` | `plink_pattern` | `data.plink_pattern` | the default is the pattern of the shipped example, `"chr{chr}"` |
| `data_prep_base.R` | `gds_pattern` | `data_prep.gds_pattern` |  |
| `data_prep_base.R` | `favor_gds_pattern` | `data_prep.favor_gds_pattern` |  |
| `data_prep_base.R` | `favor_csv_pattern` | `data_prep.favor_csv_pattern` |  |
| `data_prep_base.R` | `plink_to_gds_opts` | `data_prep.plink_to_gds_opts` |  |
| `data_prep_base.R` | `favor_match_method` | `shared.favor_match_method` | one setting for the steps that share it |
| `data_prep_base.R` | `favor_features` | `data_prep.favor_features` |  |
| `data_prep_base.R` | `favor_use_xsv` | `shared.favor_use_xsv` | one setting for the steps that share it |
| `data_prep_base.R` | `favor_na_handling` | `shared.favor_na_handling` | one setting for the steps that share it |
| `data_prep_base.R` | `favor_db_format` | `shared.favor_db_format` | one setting for the steps that share it |
| `data_prep_base.R` | `favor_release` | `shared.favor_release` | one setting for the steps that share it |
| `data_prep_base.R` | `favor_rsid_policy` | `shared.favor_rsid_policy` | one setting for the steps that share it |
| `data_prep_base.R` | `favor_input_gds_dir` | `data_prep.favor_input_gds_dir` |  |
| `data_prep_base.R` | `favor_output_dir` | none | gone: the annotated files are written to `data/gds_favor/<match>/` of the dataset |
| `data_prep_base.R` | `pc_source` | `data_prep.pc_source` |  |
| `data_prep_base.R` | `n_pcs` | `data_prep.n_pcs` |  |
| `data_prep_base.R` | `pc_maf_threshold` | `data_prep.pc_maf_threshold` |  |
| `data_prep_base.R` | `pc_missing_rate` | `data_prep.pc_missing_rate` |  |
| `data_prep_base.R` | `pc_ld_threshold` | `data_prep.pc_ld_threshold` |  |
| `data_prep_base.R` | `pc_num_thread` | `data_prep.pc_num_thread` |  |
| `data_prep_base.R` | `pc_seed` | `data_prep.pc_seed` |  |
| `data_prep_base.R` | `trait` | `data_prep.trait` |  |
| `data_prep_base.R` | `sample_gds_chr` | `data_prep.sample_gds_chr` |  |
| `data_prep_base.R` | `drop_incomplete` | `data_prep.drop_incomplete` |  |
| `data_prep_base.R` | `pheno_out_name` | `data_prep.pheno_out_name` |  |
| `training_base.R` | `symlinked_shared_root` | none | gone: `output_root` says where the outputs go |
| `training_base.R` | `b_trait_type` | `training.b_trait_type` |  |
| `training_base.R` | `b_format` | `training.b_format` |  |
| `training_base.R` | `b_column_mapping` | `training.b_column_mapping` |  |
| `training_base.R` | `b_filter_include` | `training.b_filter_include` |  |
| `training_base.R` | `b_filter_exclude` | `training.b_filter_exclude` |  |
| `training_base.R` | `b_filter_pattern` | `training.b_filter_pattern` |  |
| `training_base.R` | `b_filter_autosomes` | `training.b_filter_autosomes` |  |
| `training_base.R` | `b_qc_filters` | `training.b_qc_filters` |  |
| `training_base.R` | `b_training_rds` | `training.b_training_rds` |  |
| `training_base.R` | `b_training_trait` | `training.b_training_trait` |  |
| `training_base.R` | `b_method` | `training.b_method` |  |
| `training_base.R` | `b_maf_col` | `training.b_maf_col` |  |
| `training_base.R` | `b_beta_col` | `training.b_beta_col` |  |
| `training_base.R` | `b_p_col` | `training.b_p_col` |  |
| `training_base.R` | `b_n_col` | `training.b_n_col` |  |
| `training_base.R` | `b_target_maf` | `training.b_target_maf` |  |
| `training_base.R` | `b_target_trait` | `training.b_target_trait` |  |
| `training_base.R` | `b_target_case_prop` | `training.b_target_case_prop` |  |
| `training_base.R` | `pi_format` | `training.pi_format` |  |
| `training_base.R` | `pi_column_mapping` | `training.pi_column_mapping` |  |
| `training_base.R` | `pi_filter_include` | `training.pi_filter_include` |  |
| `training_base.R` | `pi_filter_exclude` | `training.pi_filter_exclude` |  |
| `training_base.R` | `pi_exclude_authors` | `training.pi_exclude_authors` |  |
| `training_base.R` | `pi_exclude_chromosomes` | `training.pi_exclude_chromosomes` |  |
| `training_base.R` | `pi_qc_filters` | `training.pi_qc_filters` |  |
| `training_base.R` | `pi_case_csv` | `training.pi_case_csv` |  |
| `training_base.R` | `favor_match_method` | `shared.favor_match_method` | one setting for the steps that share it |
| `training_base.R` | `favor_features` | `training.favor_features` |  |
| `training_base.R` | `favor_use_xsv` | `shared.favor_use_xsv` | one setting for the steps that share it |
| `training_base.R` | `favor_na_handling` | `shared.favor_na_handling` | one setting for the steps that share it |
| `training_base.R` | `favor_db_format` | `shared.favor_db_format` | one setting for the steps that share it |
| `training_base.R` | `favor_release` | `shared.favor_release` | one setting for the steps that share it |
| `training_base.R` | `favor_rsid_policy` | `shared.favor_rsid_policy` | one setting for the steps that share it |
| `training_base.R` | `pi_case_annotated_csv` | `training.pi_case_annotated_csv` |  |
| `training_base.R` | `pi_n_models` | `training.pi_n_models` |  |
| `training_base.R` | `pi_controls_multiplier` | `training.pi_controls_multiplier` |  |
| `training_base.R` | `pi_max_controls` | `training.pi_max_controls` |  |
| `training_base.R` | `pi_model_type` | `training.pi_model_type` |  |
| `training_base.R` | `pi_chromosomes` | `training.pi_chromosomes` |  |
| `training_base.R` | `pi_features` | `shared.pi_features` | one setting for the steps that share it |
| `training_base.R` | `pi_random_seed` | `training.pi_random_seed` |  |
| `training_base.R` | `pi_models_subdir` | `training.pi_models_subdir` |  |
| `training_base.R` | `pi_models_dir` | `training.pi_models_dir` |  |
| `training_base.R` | `pi_control_format` | `training.pi_control_format` |  |
| `training_base.R` | `pi_eval_max_controls` | `training.pi_eval_max_controls` |  |
| `training_base.R` | `pi_n_eval_controls` | `training.pi_n_eval_controls` |  |
| `single_variant_base.R` | `symlinked_shared_root` | none | gone: `output_root` says where the outputs go |
| `single_variant_base.R` | `chroms` | `shared.chroms` | one setting for the steps that share it |
| `single_variant_base.R` | `use_SPA` | `single_variant.use_spa` | renamed from `use_SPA` |
| `single_variant_base.R` | `chunk_size` | `single_variant.chunk_size` |  |
| `single_variant_base.R` | `mac_cutoff` | `single_variant.mac_cutoff` |  |
| `single_variant_base.R` | `missing_imputation` | `single_variant.missing_imputation` |  |
| `single_variant_base.R` | `n_cores` | `shared.n_cores` | one setting for the steps that share it |
| `single_variant_base.R` | `ld_window` | `single_variant.ld_window` |  |
| `single_variant_base.R` | `ld_segment` | `single_variant.ld_segment` |  |
| `single_variant_base.R` | `ld_sample_n` | `single_variant.ld_sample_n` |  |
| `single_variant_base.R` | `ld_ref_maf_min` | `single_variant.ld_ref_maf_min` |  |
| `single_variant_base.R` | `calibration_method` | `single_variant.calibration_method` |  |
| `single_variant_base.R` | `lambda_mac_min` | `single_variant.lambda_mac_min` |  |
| `snv_set_base.R` | `region_type` | `snv_set.region_type` |  |
| `snv_set_base.R` | `chroms` | `shared.chroms` | one setting for the steps that share it |
| `snv_set_base.R` | `filter_spec` | `snv_set.filter_spec` | the four fields of the base file keep their values; three more fields of `GLOWr::variant_filter()` can be set, with its defaults |
| `snv_set_base.R` | `ld_threshold` | `snv_set.ld_threshold` |  |
| `snv_set_base.R` | `mac_threshold` | `snv_set.mac_threshold` |  |
| `snv_set_base.R` | `collapse_method` | `snv_set.collapse_method` |  |
| `snv_set_base.R` | `n_cores` | `shared.n_cores` | one setting for the steps that share it |
| `snv_set_base.R` | `pi_features` | `shared.pi_features` | one setting for the steps that share it |
| `snv_set_base.R` | `use_spa` | `snv_set.use_spa` |  |
| `snv_set_base.R` | `debug_errors` | `snv_set.debug_errors` |  |
| `snv_set_base.R` | `symlinked_shared_root` | none | gone: `output_root` says where the outputs go |
| `snv_set_base.R` | `calibration` | `snv_set.calibration` |  |
| `snv_set_base.R` | `staar_enabled` | `snv_set.staar_enabled` |  |
| `snv_set_base.R` | `staar_anno_features` | `snv_set.staar_anno_features` |  |
| `snv_set_base.R` | `gene_max_records` | `snv_set.gene_max_records` |  |
| `snv_set_base.R` | `gene_segment_bp` | `snv_set.gene_segment_bp` |  |
| `snv_set_base.R` | `window_size` | `snv_set.window_size` |  |
| `snv_set_base.R` | `step_size` | `snv_set.step_size` |  |
| `snv_set_base.R` | `chunk_strategy` | `snv_set.chunk_strategy` |  |
| `snv_set_base.R` | `chunk_size_mb` | `snv_set.chunk_size_mb` |  |
| `snv_set_base.R` | `staar_native` | `snv_set.staar_native` |  |
| `snv_set_base.R` | `staar_native_use_annotation` | `snv_set.staar_native_use_annotation` |  |
| `snv_set_base.R` | `categories` | `snv_set.categories` |  |
| `snv_set_base.R` | `staar_modes` | `snv_set.staar_modes` |  |
| `snv_set_base.R` | `staar_bakein_glow` | `snv_set.staar_bakein_glow` |  |
| `snv_set_base.R` | `gene_num_in_array` | `snv_set.gene_num_in_array` |  |
| `snv_set_base.R` | `Annotation_dir` | `snv_set.Annotation_dir` |  |
| `snv_set_base.R` | `staar_rv_num_cutoff` | `snv_set.staar_rv_num_cutoff` |  |
| `snv_set_base.R` | `geno_missing_imputation` | `snv_set.geno_missing_imputation` |  |
| `snv_set_base.R` | `staar_write_subtests` | `snv_set.staar_write_subtests` |  |
| `snv_set_base.R` | `noncoding_categories` | `snv_set.noncoding_categories` |  |
| `snv_set_base.R` | `write_flat_table` | `snv_set.write_flat_table` |  |
| `snv_set_base.R` | `write_evidence` | `snv_set.write_evidence` |  |
| `snv_set_base.R` | `write_staar_detail` | `snv_set.write_staar_detail` |  |
| `snv_set_base.R` | `write_result_objects` | `snv_set.write_result_objects` |  |
| `snv_set_base.R` | `output_format` | `snv_set.output_format` |  |
| `snv_set_base.R` | `primary_test` | `snv_set.primary_test` |  |
| `snv_set_base.R` | `primary_test_glow` | `snv_set.primary_test_glow` |  |
| `snv_set_base.R` | `alpha` | `snv_set.alpha` |  |
| `snv_set_base.R` | `merge_gap` | `snv_set.merge_gap` |  |
| `snv_set_base.R` | `cmac_cutoff` | `snv_set.cmac_cutoff` |  |
| `snv_set_base.R` | `top_k` | `snv_set.top_k` |  |
| `snv_set_base.R` | `staar_genomewide_alpha` | `snv_set.staar_genomewide_alpha` |  |
| `data_prep_base.R` | `data_root` | none | replaced by `dataset`, `output_root` and the paths of the `data` section |
| `data_prep_base.R` | `base_name` | none | replaced by `dataset` |
| `data_prep_base.R` | `favor_db` | `data.favor_db` | one setting for the steps that share it; was unset; the default is the file of the shipped example |
| `data_prep_base.R` | `pheno_csv` | `data.pheno_csv` | was unset; the default is the file of the shipped example |
| `data_prep_base.R` | `pheno_id_col` | `data.pheno_id_col` | was unset; the default is the value of the shipped example |
| `data_prep_base.R` | `outcome` | `data_prep.outcome` | was unset; the default is the value of the shipped example |
| `data_prep_base.R` | `covariates` | `data_prep.covariates` | was unset; the default is the value of the shipped example |
| `data_prep_base.R` | `pcs` | `data_prep.pcs` | was unset; a null `path` takes the table of the PC part, a null `cols` takes the first `n_pcs` components, and `id_col` is the example's |
| `training_base.R` | `b_raw_path` | `training.b_raw_path` | was unset; the default is the file of the shipped example |
| `training_base.R` | `pi_raw_path` | `training.pi_raw_path` | was unset; the default is the file of the shipped example |
| `training_base.R` | `favor_db` | `data.favor_db` | one setting for the steps that share it; was unset; the default is the file of the shipped example |
| `training_base.R` | `pi_control_source` | `training.pi_control_source` | was unset; when null, the output of the earlier step is used |
| `single_variant_base.R` | `pheno_path` | `shared.pheno_path` | one setting for the steps that share it; was unset; when null, the output of the earlier step is used |
| `single_variant_base.R` | `gds_dir` | `shared.analysis_gds_dir` | renamed from `gds_dir`; one setting for the steps that share it; was unset; when null, the output of the earlier step is used |
| `single_variant_base.R` | `gds_pattern` | `shared.analysis_gds_pattern` | renamed from `gds_pattern`; one setting for the steps that share it; was unset; when null, the output of the earlier step is used |
| `single_variant_base.R` | `ld_scores_path` | `single_variant.ld_scores_path` | was unset; when null, the output of the earlier step is used |
| `snv_set_base.R` | `gds_dir` | `shared.analysis_gds_dir` | renamed from `gds_dir`; one setting for the steps that share it; was unset; when null, the output of the earlier step is used |
| `snv_set_base.R` | `gds_pattern` | `shared.analysis_gds_pattern` | renamed from `gds_pattern`; one setting for the steps that share it; was unset; when null, the output of the earlier step is used |
| `snv_set_base.R` | `pheno_path` | `shared.pheno_path` | one setting for the steps that share it; was unset; when null, the output of the earlier step is used |
| `snv_set_base.R` | `pi_model_dir` | `snv_set.pi_model_dir` | was unset; when null, the output of the earlier step is used |
| `snv_set_base.R` | `b_func` | `snv_set.b_func` | was unset; null by default |
| `snv_set_base.R` | `b_model_path` | `snv_set.b_model_path` | was unset; when null, the output of the earlier step is used |
| (no base file) | `output_dir` | none | gone: the output folder is derived, never set |
| (no base file) | `inflation_factor` | `snv_set.inflation_factor` | was read with an in-code default of NULL; a registered setting now, null by default |

## 11. Terra

This section is for a user of the Terra workflow.

Nothing is migrated on Terra. The inputs of `cloud/terra/glow_scan.wdl`, their names and
their types are those of 0.3.6, and so are their defaults apart from the default image. An
inputs file written for 0.3.6 runs unchanged, with one exception: a line that sets
`GLOWScan.docker_image` to an image of 0.3.6 or earlier is removed, or changed to the image of
this version. Two things differ once the new workflow version is selected. Two delivered files
have new names, `single_variant_all.csv` and `single_variant_plots.pdf`, while the workflow
outputs that deliver them keep theirs. And the `shared_artifacts` archive holds the YAML
configuration of the tasks where it held `config.R`. A run on the workflow version of 0.3.6 is
not touched. `cloud/README.md` has the details.

