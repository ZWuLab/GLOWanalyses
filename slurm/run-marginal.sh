#!/bin/bash
# ============================================================================
# slurm/run-marginal.sh - SLURM array template for the single-variant scan (M1)
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG). The
# R driver derives output_dir from it, so one template serves every run. One task
# per chromosome (--array=1-22); each task scans that chromosome
# (02-single-variant/marginal-scan.R --chr) -> results/marginal_chr<N>.csv. After
# the array finishes, run `marginal-scan.R --config <cfg> --combine` (a single
# job, not arrayed) to rbind the per-chr CSVs into marginal_all.csv + plots.
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). SLURM logs:
# there is deliberately NO `#SBATCH --output` (a static directive cannot point
# into a per-run outputs/slurm-logs/); pass --output/--error on the CLI.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD); the
#   config may be an absolute path or one relative to the GLOWanalyses dir:
#   cd /path/to/GLOWanalyses
#   CFG=/abs/path/to/runs/extended/config.R
#   LOGS=/abs/path/to/runs/extended/outputs/slurm-logs
#   sbatch --array=1-22 \
#          --output="$LOGS/marg_chr%a_%j.log" --error="$LOGS/marg_chr%a_%j.err" \
#          slurm/run-marginal.sh "$CFG"
#   # then, once the array is done:
#   Rscript 02-single-variant/marginal-scan.R --config "$CFG" --combine

#SBATCH --job-name=glow-marginal
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00          # dev-cluster default; override per HPC

set -euo pipefail

CONFIG="${1:-${GLOW_CONFIG:-}}"
if [ -z "${CONFIG}" ]; then
  echo "ERROR: pass the run config as the first arg, e.g." >&2
  echo "  sbatch --array=1-22 [flags] slurm/run-marginal.sh <run>/config.R" >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r_env

CHR=${SLURM_ARRAY_TASK_ID}
echo "=== GLOW marginal chr ${CHR} (config ${CONFIG}) starting at $(date) on $(hostname) ==="

# Stage scripts are GLOWanalyses-root-relative; the job CWD must be the GLOWanalyses
# directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "02-single-variant/marginal-scan.R" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 02-single-variant/marginal-scan.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi

Rscript 02-single-variant/marginal-scan.R \
    --chr "${CHR}" \
    --config "${CONFIG}"

echo "=== GLOW marginal chr ${CHR} finished at $(date) ==="
