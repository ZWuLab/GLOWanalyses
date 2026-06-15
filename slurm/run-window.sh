#!/bin/bash
# ============================================================================
# slurm/run-window.sh - SLURM array template for the sliding-window scan (Stage 2)
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG). The
# R driver derives output_dir from it. One task PER CHUNK (--array=1-N_CHUNKS);
# N_CHUNKS = nrow(<output_dir>/shared/01-genome_chunks.rds) from Stage 1.
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). SLURM logs:
# pass --output/--error on the CLI to route them into the run's slurm-logs/.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD):
#   cd /path/to/GLOWanalyses
#   CFG=runs/example/config.R
#   LOGS=runs/example/outputs/slurm-logs
#   # N_CHUNKS printed by prepare.R (or: Rscript -e 'nrow(readRDS(file.path(
#   #   dirname("'"$CFG"'"),"outputs","shared","01-genome_chunks.rds")))')
#   sbatch --array=1-"$N_CHUNKS" \
#          --output="$LOGS/chunk_%a_%j.log" --error="$LOGS/chunk_%a_%j.err" \
#          slurm/run-window.sh "$CFG"

#SBATCH --job-name=glow-window
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00          # dev-cluster default; override per HPC

set -euo pipefail

CONFIG="${1:-${GLOW_CONFIG:-}}"
if [ -z "${CONFIG}" ]; then
  echo "ERROR: pass the run config as the first arg, e.g." >&2
  echo "  sbatch --array=1-N [flags] slurm/run-window.sh <run>/config.R" >&2
  exit 1
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r_env

CHUNK=${SLURM_ARRAY_TASK_ID}
echo "=== GLOW window chunk ${CHUNK} (config ${CONFIG}) starting at $(date) on $(hostname) ==="

# Stage scripts are GLOWanalyses-root-relative; the job CWD must be the GLOWanalyses
# directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "03-snv-set/run-window.R" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 03-snv-set/run-window.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi

Rscript 03-snv-set/run-window.R \
    --chunk "${CHUNK}" \
    --config "${CONFIG}"

echo "=== GLOW window chunk ${CHUNK} finished at $(date) ==="
