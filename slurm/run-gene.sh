#!/bin/bash
# ============================================================================
# slurm/run-gene.sh - SLURM array template for the gene SNV-set scan (Stage 2)
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG). The
# R driver derives output_dir from it, so one template serves every run. One task
# per chromosome (--array=1-22); each task tests that chromosome's gene table.
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). R environment:
# activated inside the job by slurm/_job_lib.sh -- GLOW_CONDA_ENV (default r_env)
# or GLOW_RSCRIPT. SLURM logs: there is deliberately NO `#SBATCH --output` (a
# static directive cannot point into a per-run outputs/slurm-logs/); pass
# --output/--error on the CLI.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD, so the
#   stage script resolves on the compute node):
#   cd /path/to/GLOWanalyses
#   export GLOW_CONDA_ENV=r_env          # your conda env with the GLOW packages
#   CFG=runs/example/config.R
#   LOGS=runs/example/outputs/slurm-logs; mkdir -p "$LOGS"
#   sbatch --array=1-22 \
#          --output="$LOGS/gene_chr%a_%j.log" --error="$LOGS/gene_chr%a_%j.err" \
#          slurm/run-gene.sh "$CFG"

#SBATCH --job-name=glow-gene
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00          # dev-cluster default; override per HPC

set -euo pipefail

# Stage scripts + helpers are GLOWanalyses-root-relative; the job CWD must be the
# GLOWanalyses directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "03-snv-set/run-gene.R" ] || [ ! -f "slurm/_job_lib.sh" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 03-snv-set/run-gene.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi
source slurm/_job_lib.sh

CONFIG="${1:-${GLOW_CONFIG:-}}"
glow_require_config "${CONFIG}" "sbatch --array=1-22 [flags] slurm/run-gene.sh <run>/config.R"
CHR="${SLURM_ARRAY_TASK_ID:?submit as a SLURM array (sbatch --array=1-22 ...); SLURM_ARRAY_TASK_ID is unset}"

echo "=== GLOW gene chr ${CHR} (config ${CONFIG}) starting at $(date) on $(hostname) ==="
glow_activate_r

Rscript 03-snv-set/run-gene.R \
    --chr "${CHR}" \
    --config "${CONFIG}"

echo "=== GLOW gene chr ${CHR} finished at $(date) ==="
