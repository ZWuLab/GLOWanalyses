# Terra configuration files for `glow_scan.wdl`

One pair of files per run type, uploaded on the workflow's configuration page in Terra: the
inputs file binds the workflow's inputs to the columns of a data table, and the outputs file
writes each named output back to a column of the same row.

| Run type | Inputs file | Outputs file |
|---|---|---|
| single-variant scan with calibration | `GLOWScan_single_variant.inputs.json` | `GLOWScan_single_variant.outputs.json` |
| gene | `GLOWScan_gene.inputs.json` | `GLOWScan_gene.outputs.json` |
| window | `GLOWScan_window.inputs.json` | `GLOWScan_window.outputs.json` |
| coding | `GLOWScan_coding.inputs.json` | `GLOWScan_coding.outputs.json` |
| noncoding | `GLOWScan_noncoding.inputs.json` | `GLOWScan_noncoding.outputs.json` |

The files are in the format Terra's own "Download json" writes: a value that reads a column is
`${this.<column>}`, and a literal is a JSON value. The inputs read the columns `gds_files`,
`chromosomes`, `pheno_bundle` and, for a set-based run, `pi_models_tar`, `b_model_rds` and
`calibration_report`. The outputs write `glow_<run type>_<output>` columns, except that the
single-variant run's `calibration_report_out` writes the `calibration_report` column, which the
set-based runs read. The pairs carry no analysis settings, so a run keeps the workflow's
defaults, which are GLOW's full-spectrum design rather than STAARpipeline's rare-variant
conventions. To change a setting, add a line such as `"GLOWScan.rare_maf_cutoff": 0.01` to the
inputs file; `../../README.md` ("Defaults to check before a real run") lists the five lines that
switch a run to the STAARpipeline settings. The column table, the outputs per run type and the
order of the runs are in the same README ("Running on Terra").

Derived from the configuration files written by Alberto Brusati (Landers lab, UMass Chan Medical
School), 2026, for his Terra workflows.
