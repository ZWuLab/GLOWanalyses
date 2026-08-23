#!/bin/bash
# ============================================================================
# slurm/run-coding.sh - SLURM array template for the coding SNV-set scan (Stage 2)
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG). The
# R driver derives output_dir from it. One task per chromosome (--array=1-22);
# each task loops the configured coding categories over that chromosome's gene
# table. Probe chr22 first to calibrate --time (envelope ~2-4x the gene run).
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
#   sbatch --array=1-22 \
#          --output="$LOGS/coding_chr%a_%j.log" --error="$LOGS/coding_chr%a_%j.err" \
#          slurm/run-coding.sh "$CFG"

#SBATCH --job-name=glow-coding
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00          # dev-cluster default; override per HPC (probe chr22)

set -euo pipefail

# Stage scripts + helpers are GLOWanalyses-root-relative; the job CWD must be the
# GLOWanalyses directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "03-snv-set/run-coding.R" ] || [ ! -f "slurm/_job_lib.sh" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 03-snv-set/run-coding.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi
source slurm/_job_lib.sh

CONFIG="${1:-${GLOW_CONFIG:-}}"
glow_require_config "${CONFIG}" "sbatch --array=1-22 [flags] slurm/run-coding.sh <run>/config.R"
CHR="${SLURM_ARRAY_TASK_ID:?submit as a SLURM array (sbatch --array=1-22 ...); SLURM_ARRAY_TASK_ID is unset}"

echo "=== GLOW coding chr ${CHR} (config ${CONFIG}) starting at $(date) on $(hostname) ==="
glow_activate_r

Rscript 03-snv-set/run-coding.R \
    --chr "${CHR}" \
    --config "${CONFIG}"

echo "=== GLOW coding chr ${CHR} finished at $(date) ==="
