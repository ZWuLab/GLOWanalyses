#!/bin/bash
# ============================================================================
# slurm/run-annotate-favor.sh - SLURM array template for FAVOR annotation (00)
#
# RUN-AGNOSTIC: the data-prep run config is the FIRST positional arg (or
# $GLOW_CONFIG). One task per chromosome (--array=1-22, or the range your
# config's `chroms` covers); each task runs 00-data-prep/annotate-favor.R --chr N,
# i.e. <base_name>_gds/chr<N>... -> <base_name>_gds_favor/<match>/{gds,csv}/.
# FAVOR annotation is the heaviest 00 step (each FAVOR Essential-DB chunk is a
# 20-30 GB CSV), so at whole-genome scale run it as this array rather than the
# serial all-chromosome loop. Outputs are cohort data assets under data_root (not
# run outputs/); each task writes its own provenance/config_snapshot_chr<N>.rds.
# compute-pcs.R and assemble-pheno-covar.R follow serially (no array).
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). Measured on a
# WGS cohort against the FAVOR Essential DB (flexible match): ~13-18 min and a
# ~32 GB peak RSS per chromosome, hence --mem=40G; calibrate on one chromosome
# first. R environment: activated inside the job by slurm/_job_lib.sh --
# GLOW_CONDA_ENV (default r_env) or GLOW_RSCRIPT. SLURM logs: there is
# deliberately NO `#SBATCH --output` (a static directive cannot point into a
# per-run dir); pass --output/--error on the CLI.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD):
#   cd /path/to/GLOWanalyses
#   export GLOW_CONDA_ENV=r_env          # your conda env with the GLOW packages
#   CFG=runs/<name>/data-prep/config.R
#   LOGS=runs/<name>/data-prep/outputs/slurm-logs; mkdir -p "$LOGS"
#   sbatch --array=1-22 \
#          --output="$LOGS/favor_chr%a_%j.log" --error="$LOGS/favor_chr%a_%j.err" \
#          slurm/run-annotate-favor.sh "$CFG"

#SBATCH --job-name=glow-favor
#SBATCH --mem=40G                # dev-cluster default (~32 GB peak measured); override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00          # dev-cluster default (~15 min/chr measured); override per HPC

set -euo pipefail

# Stage scripts + helpers are GLOWanalyses-root-relative; the job CWD must be the
# GLOWanalyses directory (SLURM inherits the submit CWD). Fail fast with guidance otherwise.
if [ ! -f "00-data-prep/annotate-favor.R" ] || [ ! -f "slurm/_job_lib.sh" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 00-data-prep/annotate-favor.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi
source slurm/_job_lib.sh

CONFIG="${1:-${GLOW_CONFIG:-}}"
glow_require_config "${CONFIG}" "sbatch --array=1-22 [flags] slurm/run-annotate-favor.sh <run>/data-prep/config.R"
CHR="${SLURM_ARRAY_TASK_ID:?submit as a SLURM array (sbatch --array=1-22 ...); SLURM_ARRAY_TASK_ID is unset}"

echo "=== GLOW annotate-favor chr ${CHR} (config ${CONFIG}) starting at $(date) on $(hostname) ==="
glow_activate_r

Rscript 00-data-prep/annotate-favor.R \
    --config "${CONFIG}" \
    --chr "${CHR}"

echo "=== GLOW annotate-favor chr ${CHR} finished at $(date) ==="
