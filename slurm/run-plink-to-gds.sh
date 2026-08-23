#!/bin/bash
# ============================================================================
# slurm/run-plink-to-gds.sh - SLURM array template for PLINK -> GDS conversion (00)
#
# RUN-AGNOSTIC: the data-prep run config is the FIRST positional arg (or
# $GLOW_CONFIG). One task per chromosome (--array=1-22, or the range your
# config's `chroms` covers); each task runs 00-data-prep/plink-to-gds.R --chr N,
# i.e. <base_name>_plink/chr<N>.{bed,bim,fam} -> <base_name>_gds/chr<N>....gds.
# Outputs are cohort data assets under data_root (not run outputs/); each task
# writes its own provenance/config_snapshot_chr<N>.rds. Chain the FAVOR array
# after it with `sbatch --dependency=afterok:<jobid> ... slurm/run-annotate-favor.sh`.
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH); the defaults
# below are unmeasured dev-cluster values -- calibrate on one chromosome first.
# R environment: activated inside the job by slurm/_job_lib.sh -- GLOW_CONDA_ENV
# (default r_env) or GLOW_RSCRIPT. SLURM logs: there is deliberately NO
# `#SBATCH --output`; pass --output/--error on the CLI.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD):
#   cd /path/to/GLOWanalyses
#   export GLOW_CONDA_ENV=r_env          # your conda env with the GLOW packages
#   CFG=runs/<name>/data-prep/config.R
#   LOGS=runs/<name>/data-prep/outputs/slurm-logs; mkdir -p "$LOGS"
#   sbatch --array=1-22 \
#          --output="$LOGS/gds_chr%a_%j.log" --error="$LOGS/gds_chr%a_%j.err" \
#          slurm/run-plink-to-gds.sh "$CFG"

#SBATCH --job-name=glow-plink2gds
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00          # dev-cluster default; override per HPC

set -euo pipefail

# Stage scripts + helpers are GLOWanalyses-root-relative; the job CWD must be the
# GLOWanalyses directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "00-data-prep/plink-to-gds.R" ] || [ ! -f "slurm/_job_lib.sh" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 00-data-prep/plink-to-gds.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi
source slurm/_job_lib.sh

CONFIG="${1:-${GLOW_CONFIG:-}}"
glow_require_config "${CONFIG}" "sbatch --array=1-22 [flags] slurm/run-plink-to-gds.sh <run>/data-prep/config.R"
CHR="${SLURM_ARRAY_TASK_ID:?submit as a SLURM array (sbatch --array=1-22 ...); SLURM_ARRAY_TASK_ID is unset}"

echo "=== GLOW plink-to-gds chr ${CHR} (config ${CONFIG}) starting at $(date) on $(hostname) ==="
glow_activate_r

Rscript 00-data-prep/plink-to-gds.R \
    --config "${CONFIG}" \
    --chr "${CHR}"

echo "=== GLOW plink-to-gds chr ${CHR} finished at $(date) ==="
