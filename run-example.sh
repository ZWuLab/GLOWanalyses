#!/usr/bin/env bash
# ============================================================================
# run-example.sh - run the shipped example from start to finish
#
# The example is a small synthetic cohort (data-example/) with the configuration
# example.yml. This script runs it through the `glow` command and then checks that
# every step left a valid, non-empty result. It shows that the workflow runs on a
# fresh clone without any other data: the example ships its annotated GDS files,
# its phenotype bundle and its B and PI models, so the default run needs only the
# two analysis steps.
#
#   ./glow run single_variant     the single-variant scan, the LD scores, the calibration
#   ./glow run snv_set            the gene, window and coding scans, each with its summary
#
# The results are valid in form and without meaning as statistics: the data are
# synthetic. They are written to runs/example/baseline/.
#
# Usage (with GLOWr and GLOWpipeline installed; the script runs from any directory):
#   conda activate r_env
#   bash run-example.sh               # the two analysis steps and the check
#   REBUILD=1 bash run-example.sh     # first build the cohort, the data preparation and the training again
#   DRYRUN=1  bash run-example.sh     # show what each command would run, and run nothing
#
#   REBUILD=1   regenerates the synthetic cohort (data-example/make-example-cohort.R)
#               and runs the data preparation and the training with --force, which
#               replaces the pre-built files under runs/example/. Use it after a
#               change of the generator. data-example/build-intermediates.sh copies
#               the rebuilt files back into data-example/ for a release.
#   DRYRUN=1    passes --dry-run to every `glow run`: no unit is started. The pre-built
#               files of the example are still put into runs/example/, so that the
#               states shown are those of a real run.
#
# A second run skips what is finished: every unit of work keeps a completion record.
# Exit status: non-zero when a step or the check fails.
set -euo pipefail
cd "$(dirname "$0")"

if [[ ! -d data-example || ! -x glow ]]; then
  echo "ERROR: expected the GLOWanalyses directory (data-example/ and the glow command)." >&2
  exit 1
fi

REBUILD="${REBUILD:-0}"
DRYRUN="${DRYRUN:-0}"
DRY=(); if [[ "$DRYRUN" == "1" ]]; then DRY=(--dry-run); fi
say() { echo "+ $*"; }
run() { say "$@"; "$@"; }

echo "=================================================================="
echo "  GLOWanalyses example  (REBUILD=$REBUILD  DRYRUN=$DRYRUN)"
echo "=================================================================="

FORCE=()
if [[ "$REBUILD" == "1" ]]; then
  echo "--- [REBUILD] the synthetic cohort, the data preparation and the training ---"
  if [[ "$DRYRUN" == "1" ]]; then say Rscript data-example/make-example-cohort.R
  else run Rscript data-example/make-example-cohort.R; fi
  run ./glow run data_prep --force --wait "${DRY[@]+"${DRY[@]}"}"
  run ./glow run training --force --wait "${DRY[@]+"${DRY[@]}"}"
  # the check below then requires the rebuilt outputs instead of skipping them
  export GLOW_EXAMPLE_REBUILT=1
  # the steps below read what was just replaced, so their finished units are out of date
  # and run again by themselves; --force also covers a change of the generator's settings
  FORCE=(--force)
fi

# --wait matters when example.yml is changed to run.mode slurm: the command then stays until the
# jobs of the step have finished, and stops when one failed. In the local mode it changes nothing.
echo "--- the single-variant step ---"
run ./glow run single_variant --wait "${FORCE[@]+"${FORCE[@]}"}" "${DRY[@]+"${DRY[@]}"}"

echo "--- the SNV-set step: gene, window and coding ---"
run ./glow run snv_set --wait "${FORCE[@]+"${FORCE[@]}"}" "${DRY[@]+"${DRY[@]}"}"

if [[ "$DRYRUN" == "1" ]]; then
  echo "=================================================================="
  echo "  DRYRUN complete - no unit was started."
  echo "=================================================================="
  exit 0
fi

echo "--- the check of the results ---"
run Rscript data-example/check-example.R

echo "=================================================================="
echo "  Example run complete. Results: runs/example/baseline/"
echo "=================================================================="
