#!/bin/bash
# ============================================================================
# slurm/run-staar-window.sh - SLURM array template for the NATIVE STAARpipeline
# sliding-window comparison (03-snv-set/run-staar-window.R), one task per chunk
#
# RUN-AGNOSTIC: the run config is the FIRST positional arg (or $GLOW_CONFIG);
# the SPA mode (on|off) is the SECOND positional arg (or $GLOW_SPA_MODE).
# Submit ONCE PER MODE. One task PER CHUNK (--array=1-N_CHUNKS); N_CHUNKS =
# nrow(<output_dir>/shared/01-genome_chunks.rds) from Stage 1 (prepare.R with
# staar_native = TRUE).
#
# Resources are overridable on the sbatch CLI (CLI beats #SBATCH). R
# environment: activated inside the job by slurm/_job_lib.sh -- GLOW_CONDA_ENV
# (default r_env) or GLOW_RSCRIPT. SLURM logs: pass --output/--error on the CLI
# to route them into the run's slurm-logs/.
#
#   SUBMIT FROM the GLOWanalyses directory (SLURM inherits the submit CWD):
#   cd /path/to/GLOWanalyses
#   export GLOW_CONDA_ENV=r_env
#   CFG=runs/example/config.R
#   LOGS=runs/example/outputs/slurm-logs; mkdir -p "$LOGS"
#   sbatch --array=1-"$N_CHUNKS" \
#          --output="$LOGS/staar_win_off_%a_%j.log" --error="$LOGS/staar_win_off_%a_%j.err" \
#          slurm/run-staar-window.sh "$CFG" off     # non-SPA: STAAR-O + grid
#   sbatch --array=1-"$N_CHUNKS" \
#          --output="$LOGS/staar_win_on_%a_%j.log" --error="$LOGS/staar_win_on_%a_%j.err" \
#          slurm/run-staar-window.sh "$CFG" on      # SPA: STAAR-B

#SBATCH --job-name=staar-window
#SBATCH --mem=16G                # dev-cluster default; override per HPC
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00          # dev-cluster default; override per HPC

set -euo pipefail

if [ ! -f "03-snv-set/run-staar-window.R" ] || [ ! -f "slurm/_job_lib.sh" ]; then
  echo "ERROR: not in the GLOWanalyses directory (no 03-snv-set/run-staar-window.R in $PWD)." >&2
  echo "       Submit from the GLOWanalyses dir, or pass sbatch --chdir=/path/to/GLOWanalyses." >&2
  exit 1
fi
source slurm/_job_lib.sh

CONFIG="${1:-${GLOW_CONFIG:-}}"
glow_require_config "${CONFIG}" "sbatch --array=1-N [flags] slurm/run-staar-window.sh <run>/config.R {on|off}"
SPA_MODE="${2:-${GLOW_SPA_MODE:-}}"
if [ "${SPA_MODE}" != "on" ] && [ "${SPA_MODE}" != "off" ]; then
  echo "ERROR: SPA mode must be 'on' or 'off' (2nd arg or GLOW_SPA_MODE; got '${SPA_MODE:-<unset>}')." >&2
  exit 1
fi
CHUNK="${SLURM_ARRAY_TASK_ID:?submit as a SLURM array (sbatch --array=1-N ...); SLURM_ARRAY_TASK_ID is unset}"

echo "=== native STAAR window chunk ${CHUNK} (config ${CONFIG}, spa ${SPA_MODE}) starting at $(date) on $(hostname) ==="
glow_activate_r

Rscript 03-snv-set/run-staar-window.R \
    --chunk "${CHUNK}" \
    --config "${CONFIG}" \
    --spa "${SPA_MODE}"

echo "=== native STAAR window chunk ${CHUNK} (spa ${SPA_MODE}) finished at $(date) ==="
