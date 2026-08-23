#!/bin/bash
# ============================================================================
# slurm/run-window.sh - SLURM array template for the sliding-window scan (Stage 2)
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG). The
# R driver derives output_dir from it. One task PER CHUNK (--array=1-N_CHUNKS);
# N_CHUNKS = nrow(<output_dir>/shared/01-genome_chunks.rds) from Stage 1.
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). R environment:
# activated inside the job by slurm/_job_lib.sh -- GLOW_CONDA_ENV (default r_env)
# or GLOW_RSCRIPT. SLURM logs: pass --output/--error on the CLI to route them
# into the run's slurm-logs/.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD):
#   cd /path/to/GLOWanalyses
#   export GLOW_CONDA_ENV=r_env          # your conda env with the GLOW packages
#   CFG=runs/example/config.R
#   LOGS=runs/example/outputs/slurm-logs; mkdir -p "$LOGS"
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

# Stage scripts + helpers are GLOWanalyses-root-relative; the job CWD must be the
# GLOWanalyses directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "03-snv-set/run-window.R" ] || [ ! -f "slurm/_job_lib.sh" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 03-snv-set/run-window.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi
source slurm/_job_lib.sh

CONFIG="${1:-${GLOW_CONFIG:-}}"
glow_require_config "${CONFIG}" "sbatch --array=1-N [flags] slurm/run-window.sh <run>/config.R"
CHUNK="${SLURM_ARRAY_TASK_ID:?submit as a SLURM array (sbatch --array=1-N ...); SLURM_ARRAY_TASK_ID is unset}"

echo "=== GLOW window chunk ${CHUNK} (config ${CONFIG}) starting at $(date) on $(hostname) ==="
glow_activate_r

Rscript 03-snv-set/run-window.R \
    --chunk "${CHUNK}" \
    --config "${CONFIG}"

echo "=== GLOW window chunk ${CHUNK} finished at $(date) ==="
