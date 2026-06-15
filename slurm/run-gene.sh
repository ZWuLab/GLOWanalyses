#!/bin/bash
# ============================================================================
# slurm/run-gene.sh - SLURM array template for the gene SNV-set scan (Stage 2)
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG). The
# R driver derives output_dir from it, so one template serves every run. One task
# per chromosome (--array=1-22); each task tests that chromosome's gene table.
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). SLURM logs:
# there is deliberately NO `#SBATCH --output` (a static directive cannot point
# into a per-run outputs/slurm-logs/); pass --output/--error on the CLI.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD, so the
#   stage script resolves on the compute node):
#   cd /path/to/GLOWanalyses
#   CFG=runs/example/config.R
#   LOGS=runs/example/outputs/slurm-logs
#   sbatch --array=1-22 \
#          --output="$LOGS/gene_chr%a_%j.log" --error="$LOGS/gene_chr%a_%j.err" \
#          slurm/run-gene.sh "$CFG"

#SBATCH --job-name=glow-gene
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00          # dev-cluster default; override per HPC

set -euo pipefail

CONFIG="${1:-${GLOW_CONFIG:-}}"
if [ -z "${CONFIG}" ]; then
  echo "ERROR: pass the run config as the first arg, e.g." >&2
  echo "  sbatch --array=1-22 [flags] slurm/run-gene.sh <run>/config.R" >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r_env

CHR=${SLURM_ARRAY_TASK_ID}
echo "=== GLOW gene chr ${CHR} (config ${CONFIG}) starting at $(date) on $(hostname) ==="

# Stage scripts are GLOWanalyses-root-relative; the job CWD must be the GLOWanalyses
# directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "03-snv-set/run-gene.R" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 03-snv-set/run-gene.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi

Rscript 03-snv-set/run-gene.R \
    --chr "${CHR}" \
    --config "${CONFIG}"

echo "=== GLOW gene chr ${CHR} finished at $(date) ==="
