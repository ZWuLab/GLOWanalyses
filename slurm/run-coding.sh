#!/bin/bash
# ============================================================================
# slurm/run-coding.sh - SLURM array template for the coding SNV-set scan (Stage 2)
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG). The
# R driver derives output_dir from it. One task per chromosome (--array=1-22);
# each task loops the configured coding categories over that chromosome's gene
# table. Probe chr22 first to calibrate --time (envelope ~2-4x the gene run).
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). SLURM logs:
# pass --output/--error on the CLI to route them into the run's slurm-logs/.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD):
#   cd /path/to/GLOWanalyses
#   CFG=runs/example/config.R
#   LOGS=runs/example/outputs/slurm-logs
#   sbatch --array=1-22 \
#          --output="$LOGS/coding_chr%a_%j.log" --error="$LOGS/coding_chr%a_%j.err" \
#          slurm/run-coding.sh "$CFG"

#SBATCH --job-name=glow-coding
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00          # dev-cluster default; override per HPC (probe chr22)

set -euo pipefail

CONFIG="${1:-${GLOW_CONFIG:-}}"
if [ -z "${CONFIG}" ]; then
  echo "ERROR: pass the run config as the first arg, e.g." >&2
  echo "  sbatch --array=1-22 [flags] slurm/run-coding.sh <run>/config.R" >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r_env

CHR=${SLURM_ARRAY_TASK_ID}
echo "=== GLOW coding chr ${CHR} (config ${CONFIG}) starting at $(date) on $(hostname) ==="

# Stage scripts are GLOWanalyses-root-relative; the job CWD must be the GLOWanalyses
# directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "03-snv-set/run-coding.R" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 03-snv-set/run-coding.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi

Rscript 03-snv-set/run-coding.R \
    --chr "${CHR}" \
    --config "${CONFIG}"

echo "=== GLOW coding chr ${CHR} finished at $(date) ==="
