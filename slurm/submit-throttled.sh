#!/usr/bin/env bash
# ============================================================================
# slurm/submit-throttled.sh — drip-feed a Stage-2 array under a SLURM submit cap.
#
# RUN-AGNOSTIC + REGION-TYPE-AWARE. The run config is the FIRST positional arg.
# The submitter reads `region_type` from the config and picks the matching array
# script + array dimension automatically:
#   gene   -> slurm/run-gene.sh    , one task per chromosome (default 1..max(chroms) or 22)
#   coding -> slurm/run-coding.sh  , one task per chromosome (default 1..max(chroms) or 22)
#   window -> slurm/run-window.sh  , one task per chunk (default 1..nrow(01-genome_chunks.rds))
# It derives the run's slurm-logs/ from the config and passes the config (and
# --output/--error) to the array script, so logs land in the RUN's
# outputs/slurm-logs/.
#
# WHY THIS EXISTS
#   Many clusters cap how many jobs a user may have queued at once via the QOS
#   limit `MaxSubmitJobsPerUser`. A SLURM *array* counts each task individually
#   toward that cap, so one big `sbatch --array=1-565 …` is rejected with
#   `QOSMaxSubmitJobPerUserLimit`. The `%N` array throttle does NOT help — it
#   limits how many tasks RUN at once, but all elements still count as
#   *submitted*, so the submit-time check still rejects them. (Find your cap with
#     sacctmgr -n show qos format=Name,MaxSubmitJobsPerUser )
#
# WHAT IT DOES
#   Submits units [FIRST..LAST] in waves: each pass counts the user's active jobs
#   (`squeue`), computes the room left under CAP, submits one `--array=<a>-<b>`
#   batch sized to fill that room, sleeps POLL seconds, and repeats — keeping
#   ~CAP tasks in flight without exceeding the cap. It THROTTLES SUBMISSION ONLY;
#   verify completion afterward via the per-unit logs (<run>/outputs/logs/),
#   the SLURM logs (<run>/outputs/slurm-logs/), or `sacct`, and re-submit any
#   failed unit ids.
#
# REQUIREMENTS
#   - `sbatch` + `squeue` on PATH (a SLURM submit host). No conda/R needed for
#     submission itself; the array script activates r_env inside each job.
#   - `Rscript` is used ONLY to read `region_type`/`chroms` from the config and to
#     derive a default LAST. Pass REGION_TYPE and LAST explicitly to avoid R.
#   - Run it so it survives logout — under tmux/screen, or with nohup:
#       nohup bash slurm/submit-throttled.sh \
#             runs/example/config.R > submit.log 2>&1 &
#
# USAGE (from the GLOWanalyses directory; the submitter cd's to it via $0)
#   CFG=runs/example/config.R
#   # whole run (FIRST defaults to 1; LAST derived from region_type):
#   bash slurm/submit-throttled.sh "$CFG"
#   # one unit, to calibrate wall-time:
#   bash slurm/submit-throttled.sh "$CFG" 1 1 --time=02:00:00
#   # a unit-id range:
#   bash slurm/submit-throttled.sh "$CFG" 1 22
#   # tune the cap + poll, preview without submitting:
#   CAP=90 POLL=180 bash .../submit-throttled.sh "$CFG"
#   DRYRUN=1 bash .../submit-throttled.sh "$CFG"
#
# CONFIG (environment variables)
#   CAP          max jobs to keep in the queue (default 90; stay < your cap)
#   POLL         seconds between queue checks (default 120)
#   DRYRUN       if set, print the sbatch command for every batch it WOULD submit
#                (assuming an empty queue) and exit without submitting or sleeping
#   REGION_TYPE  gene|window|coding — overrides reading it from the config
#   SCRIPT       the Stage-2 array script — overrides the region_type default
#   CHUNK_TABLE  window only; Stage-1 chunk table used to derive LAST
#                (default <run>/outputs/shared/01-genome_chunks.rds)
#   GLOW_CONFIG  alternative to the positional config arg
#   Extra sbatch flags pass through, e.g.:  CAP=90 bash .../submit-throttled.sh "$CFG" 1 22 --mem=32G

set -uo pipefail   # NOT -e: the loop must tolerate a transient sbatch/squeue miss

# Anchor to the GLOWanalyses root (this submitter lives in slurm/) so the array
# SCRIPT path + any GLOWanalyses-relative CONFIG resolve, and the array jobs inherit
# the GLOWanalyses dir as their CWD (the stage scripts are root-relative). Runs on
# the submit host, where $0 is valid.
cd "$(dirname "$0")/.."

CAP="${CAP:-90}"
POLL="${POLL:-120}"
DRYRUN="${DRYRUN:-}"

# Positional: <config> [FIRST] [LAST] [extra sbatch flags...]
CONFIG="${1:-${GLOW_CONFIG:-}}"
if [ -z "${CONFIG}" ]; then
  echo "ERROR: first arg must be the run config, e.g." >&2
  echo "  bash $0 runs/example/config.R" >&2
  exit 1
fi
[ ! -f "${CONFIG}" ] && { echo "ERROR: run config not found: ${CONFIG}" >&2; exit 1; }
shift 1

# Run paths (output_dir is structural: dirname(config)/outputs).
RUN_DIR="$(dirname "${CONFIG}")"
OUTPUT_DIR="${RUN_DIR}/outputs"
SLURM_LOGS="${OUTPUT_DIR}/slurm-logs"

# region_type: env override wins; else read from the config via R.
REGION_TYPE="${REGION_TYPE:-}"
if [ -z "${REGION_TYPE}" ] && command -v Rscript >/dev/null 2>&1; then
  REGION_TYPE="$(Rscript -e "e<-new.env(); source('${CONFIG}', local=e); cat(if (!is.null(e\$region_type)) e\$region_type else '')" 2>/dev/null)"
fi
if [ -z "${REGION_TYPE}" ]; then
  echo "ERROR: could not determine region_type from ${CONFIG}." >&2
  echo "       Set it explicitly: REGION_TYPE=gene|window|coding bash $0 <config>" >&2
  exit 1
fi

# Array script + per-task log prefix by region_type.
case "${REGION_TYPE}" in
  gene)   DEF_SCRIPT="slurm/run-gene.sh";   PREFIX="gene_chr"   ;;
  coding) DEF_SCRIPT="slurm/run-coding.sh"; PREFIX="coding_chr" ;;
  window) DEF_SCRIPT="slurm/run-window.sh"; PREFIX="chunk"      ;;
  *) echo "ERROR: unknown region_type '${REGION_TYPE}' (expected gene|window|coding)." >&2; exit 1 ;;
esac
SCRIPT="${SCRIPT:-${DEF_SCRIPT}}"

FIRST="${1:-1}"
LAST="${2:-}"
shift "$(( $# >= 2 ? 2 : $# ))"   # remaining args (if any) pass through to sbatch
EXTRA_SBATCH=("$@")

# Derive LAST by region_type when not given.
if [ -z "${LAST}" ]; then
  if [ "${REGION_TYPE}" = "window" ]; then
    CHUNK_TABLE="${CHUNK_TABLE:-${OUTPUT_DIR}/shared/01-genome_chunks.rds}"
    if command -v Rscript >/dev/null 2>&1 && [ -f "${CHUNK_TABLE}" ]; then
      LAST="$(Rscript -e "cat(nrow(readRDS('${CHUNK_TABLE}')))" 2>/dev/null)"
    fi
    if ! [[ "${LAST}" =~ ^[0-9]+$ ]]; then
      echo "ERROR: could not derive LAST (n chunks) from ${CHUNK_TABLE}." >&2
      echo "       Run Stage 1 (prepare.R) first, or pass LAST explicitly." >&2
      exit 1
    fi
  else
    # gene/coding: array task id == chromosome; default to max(chroms) or 22.
    if command -v Rscript >/dev/null 2>&1; then
      LAST="$(Rscript -e "e<-new.env(); source('${CONFIG}', local=e); cat(if (!is.null(e\$chroms)) max(as.integer(e\$chroms)) else 22L)" 2>/dev/null)"
    fi
    [ -z "${LAST}" ] && LAST=22
    if ! [[ "${LAST}" =~ ^[0-9]+$ ]]; then
      echo "ERROR: could not derive LAST (max chromosome) from ${CONFIG}; pass it explicitly." >&2
      exit 1
    fi
  fi
fi

# Validate the range + the array script.
if ! [[ "${FIRST}" =~ ^[0-9]+$ && "${LAST}" =~ ^[0-9]+$ ]] || [ "${FIRST}" -gt "${LAST}" ]; then
  echo "ERROR: need integer FIRST<=LAST (got FIRST=${FIRST}, LAST=${LAST})." >&2
  exit 1
fi
[ ! -f "${SCRIPT}" ] && { echo "ERROR: array script not found: ${SCRIPT}" >&2; exit 1; }

# sbatch needs the log dir to exist before it opens --output/--error.
mkdir -p "${SLURM_LOGS}"

echo "Throttled submit: run=${RUN_DIR} | region_type=${REGION_TYPE}"
echo "  units ${FIRST}-${LAST} | CAP=${CAP} | POLL=${POLL}s | script=${SCRIPT}"
echo "  config=${CONFIG}"
echo "  slurm logs -> ${SLURM_LOGS}/"
[ "${#EXTRA_SBATCH[@]}" -gt 0 ] && echo "  extra sbatch flags: ${EXTRA_SBATCH[*]}"

# Wave loop. DRYRUN prints the planned sbatch commands (assuming an empty queue)
# without submitting or sleeping — a safe preview / smoke test.
c="${FIRST}"
while [ "${c}" -le "${LAST}" ]; do
  if [ -n "${DRYRUN}" ]; then inq=0; else inq="$(squeue -u "${USER}" -h -r 2>/dev/null | wc -l)"; fi
  room="$(( CAP - inq ))"
  if [ "${room}" -ge 1 ]; then
    hi="$(( c + room - 1 ))"
    [ "${hi}" -gt "${LAST}" ] && hi="${LAST}"
    cmd=(sbatch --array="${c}-${hi}"
         --output="${SLURM_LOGS}/${PREFIX}%a_%j.log"
         --error="${SLURM_LOGS}/${PREFIX}%a_%j.err"
         "${EXTRA_SBATCH[@]+"${EXTRA_SBATCH[@]}"}" "${SCRIPT}" "${CONFIG}")
    if [ -n "${DRYRUN}" ]; then
      echo "[DRYRUN] ${cmd[*]}"
      c="$(( hi + 1 ))"
    elif "${cmd[@]}"; then
      echo "[$(date +%H:%M:%S)] submitted ${c}-${hi}  (queue was ${inq}, room ${room})"
      c="$(( hi + 1 ))"
    else
      echo "[$(date +%H:%M:%S)] sbatch failed for ${c}-${hi} (queue ${inq}); will retry"
    fi
  else
    echo "[$(date +%H:%M:%S)] queue at ${inq} (>= CAP ${CAP}); waiting ${POLL}s..."
  fi
  if [ "${c}" -le "${LAST}" ] && [ -z "${DRYRUN}" ]; then sleep "${POLL}"; fi
done

echo "Done: units ${FIRST}-${LAST} all submitted. Verify completion via ${OUTPUT_DIR}/logs/ and sacct."
