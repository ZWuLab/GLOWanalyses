#!/usr/bin/env bash
# ============================================================================
# run-example.sh - run the self-contained synthetic example end-to-end
#
# Drives the bundled synthetic-cohort demo (data-example/ + runs/example/) over
# the GLOWanalyses workflow templates, then runs the schema self-test. It proves
# the templates work out-of-the-box on a fresh clone with NO external data: the
# committed cohort already ships the pre-built aGDS + B/PI models, so the default
# run executes only the ANALYSIS stages 02 -> 03 -> 04 over those tracked assets.
#
# The outputs are MECHANICALLY valid but STATISTICALLY meaningless (synthetic
# data) - the example validates template MECHANICS, not statistical correctness.
#
# Usage (with GLOWr + GLOWpipeline installed; the script cd's to its own dir, so
# paths are GLOWanalyses-root-relative and it runs from anywhere):
#   conda activate r_env
#   bash run-example.sh        # 02 -> 03 -> 04 + self-test
#   REBUILD=1 bash run-example.sh   # also regen cohort + 00/01
#   DRYRUN=1  bash run-example.sh   # print the chain, run nothing
#
# Env toggles:
#   REBUILD=1   first regenerate the synthetic cohort (make-example-cohort.R) and
#               the pre-built 00/01 intermediates (build-intermediates.sh), then
#               run 02 -> 03 -> 04. Use after editing the generator.
#   DRYRUN=1    print every command without executing (a preview of the chain).
#
# Exit status: non-zero if any stage or the self-test fails.
set -euo pipefail

# Anchor to the GLOWanalyses root (this script's dir) so every stage/config path is
# GLOWanalyses-root-relative and the demo runs from any CWD. The configs source
# base-config/ + data-example/ via these same root-relative paths.
cd "$(dirname "$0")"
A="."
DE="$A/data-example"
RUN="$A/runs/example"

if [[ ! -d "$DE" ]]; then
  echo "ERROR: expected to be in the GLOWanalyses root (no $DE found)." >&2
  exit 1
fi

REBUILD="${REBUILD:-0}"
DRYRUN="${DRYRUN:-0}"

# run/say: echo each command; execute it unless DRYRUN=1.
say() { echo "+ $*"; }
run() { say "$@"; if [[ "$DRYRUN" != "1" ]]; then "$@"; fi; }

echo "=================================================================="
echo "  GLOWanalyses synthetic example  (REBUILD=$REBUILD  DRYRUN=$DRYRUN)"
echo "=================================================================="

# ---- Optional: regenerate the cohort + the pre-built 00/01 intermediates -----
if [[ "$REBUILD" == "1" ]]; then
  echo "--- [REBUILD] regenerate synthetic cohort + 00/01 intermediates ---"
  run Rscript "$DE/make-example-cohort.R"
  # build-intermediates.sh itself runs 00/plink-to-gds -> annotate-favor ->
  # compute-pcs -> assemble-pheno-covar and 01-training/{b,pi}, then copies the
  # B model + PI ensemble + pheno bundle up to the cohort root (tracked assets).
  if [[ "$DRYRUN" == "1" ]]; then
    say bash "$DE/build-intermediates.sh"
  else
    bash "$DE/build-intermediates.sh"
  fi
fi

# ---- Stage 02: single-variant (marginal-scan -> compute-ld-scores -> calibrate)
echo "--- [02] single-variant ---"
SV="$RUN/single-variant/config.R"
run Rscript "$A/02-single-variant/marginal-scan.R"     --config "$SV"
run Rscript "$A/02-single-variant/compute-ld-scores.R" --config "$SV"
run Rscript "$A/02-single-variant/calibrate.R"         --config "$SV"

# ---- Stage 03: SNV-set, gene style (prepare once, then one unit per chromosome)
echo "--- [03] SNV-set (gene) ---"
SET="$RUN/snv-set/config.R"
run Rscript "$A/03-snv-set/prepare.R" --config "$SET"
for CHR in 21 22; do
  run Rscript "$A/03-snv-set/run-gene.R" --config "$SET" --chr "$CHR"
done

# ---- Stage 04: summary (aggregate -> plots) ---------------------------------
echo "--- [04] summary ---"
run Rscript "$A/04-summary/aggregate.R" --config "$SET"
run Rscript "$A/04-summary/plots.R"     --config "$SET"

# ---- Self-test: presence + schema + non-emptiness per stage -----------------
echo "--- [check] schema self-test ---"
run Rscript "$DE/check-example.R"

if [[ "$DRYRUN" == "1" ]]; then
  echo "=================================================================="
  echo "  DRYRUN complete - no commands were executed."
  echo "=================================================================="
else
  echo "=================================================================="
  echo "  Example run complete."
  echo "=================================================================="
fi
