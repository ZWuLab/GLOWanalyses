
# `check-workflow-defaults.R`: the Terra workflow against the settings registry

## What it does

GLOWanalyses runs an analysis in two ways: the `glow` command with a YAML configuration, and
the Terra workflow `cloud/terra/glow_scan.wdl` with typed inputs. Both have to give the same
results for the same settings. They share one source of defaults, the settings registry
(`base-config/settings-registry.yml`). This script checks that the workflow does not carry a
default of its own without saying so.

It compares two things with the registry.

- **The default of every workflow input** with the default of the setting it sets.
- **Every value that a task writes into its configuration on its own**, such as the empty list
  of native comparison modes, with the default of that setting.

The file `cloud/terra/workflow-settings.yml` is the link between the two sides. It names the
setting of each input, or says why an input has none (a resource request, for example). It
also declares each deliberate difference with its reason. The script fails when

- a default differs and no reason is declared, or a reason is declared and the defaults are
  the same;
- an input is missing from the file, or the file names an input that the workflow does not
  have;
- a task writes a setting that the registry does not know, or a value that differs from the
  registry default without a declared reason;
- two input names are the same when letter case is ignored, which Terra's configuration store
  does not accept.

## How to run it

Run it from the GLOWanalyses directory, in an R environment that has GLOWpipeline with
`read_glow_config()` and the package `yaml`.

```bash
Rscript cloud/terra/check-workflow-defaults.R                       # the shipped workflow
Rscript cloud/terra/check-workflow-defaults.R --out defaults.md     # with a Markdown report
Rscript cloud/terra/check-workflow-defaults.R --wdl other.wdl       # another copy of the workflow
```

## What it prints

The report has a summary with the counts, one row per input, and one row per setting that a
task writes. The last line on the screen is `Verdict for the workflow defaults: PASS` or
`... FAIL`, and the exit status is 1 on a failure. The image build runs the script, so a
workflow default that drifts from the registry stops the build.

## When you change a default

Change it in the registry, regenerate `base-config/base_config.yml`, and change the workflow
input to the same value. When the workflow has to differ, add `differs:` with the reason to the
input's entry in `cloud/terra/workflow-settings.yml`. The reason then appears in the report.
