
# `check-glow-slurm.sh`: check the SLURM mode of the `glow` command

## What it does

With `run.mode: slurm` in the configuration, `glow run <step>` submits one job per unit of
work and returns. This script checks that mode on a cluster, with the shipped example. It
submits real jobs, about thirty in all, each of a few seconds. It makes eight groups of
checks, 35 checks in all:

- **(a)** One call of `glow run snv_set` finishes the example's gene run: the preparation, the
  two scans, the aggregation and the plots, each job waiting for the jobs it reads.
- **(b)** A second call submits nothing and reports every unit as finished.
- **(c)** After a deleted result table of one scan unit, a call submits that unit and the parts
  after it, and nothing else.
- **(d)** The window run finishes from one call with `--wait`. Its chunks are known only after
  the preparation, so the command waits for the preparation and then submits them.
- **(e)** A unit whose job fails is submitted again by the next call, together with the units
  that read it.
- **(f)** The tables of the gene and window runs equal those of a local run of the same steps,
  within a relative tolerance of 1e-8. The tolerance is there because the jobs may run on other
  processors than the local run, and the last digits of a number can then differ.
- **(h)** A call with `--wait` stops with an error when a job that it follows fails, after one
  submission of each unit, and the same command then finishes the run.
- **(g)** A job whose configuration was changed after its submission does not run. The jobs of
  this group are submitted to start 75 seconds later, the script adds a line to the
  configuration in between, and the log of the first job has to say why it stopped.

## How to run it

Run it from a shell where `sbatch` works and R is set up with GLOWr and GLOWpipeline.

```bash
bash slurm/check-glow-slurm.sh <work folder>                       # default sbatch options
bash slurm/check-glow-slurm.sh <work folder> "--partition=short"   # options for every job
```

The work folder has to be on storage that the compute nodes see, and it has to be empty or
absent. The second argument is passed to every `sbatch` call through the setting
`run.slurm.options`. The jobs set up R by activating the conda environment named in the
variable `GLOW_CONDA_ENV` (default `r_env`). The check takes five to ten minutes on an idle
cluster, and longer when the jobs wait in the queue.

## What it prints

Each check prints `[PASS]` or `[FAIL ...]` with its name. The last line is
`Verdict for the SLURM mode: PASS` or `... FAIL`, and the exit status is 1 on a failure. The
output of every command, the job scripts and the job logs are kept in the work folder.
