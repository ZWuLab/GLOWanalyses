#!/usr/bin/env bash
# ============================================================================
# cloud/terra/run-example-wdl.sh - run the synthetic example through glow_scan.wdl
#
# The WDL route's counterpart of run-example.sh: runs the bundled synthetic cohort through
# cloud/terra/glow_scan.wdl with miniwdl, seven times, then compares the tables with the ones
# run-example.sh wrote in a reference GLOWanalyses tree.
#   1. region_type = single_variant                 -> OUTDIR/single-variant
#   2. region_type = gene, no calibration           -> OUTDIR/gene             (as runs/example/snv-set)
#   3. region_type = gene, calibrated by run 1      -> OUTDIR/gene-calibrated         (single_variant_scan + ld_scores)
#   3b. the same, calibrated from run 1's report    -> OUTDIR/gene-calibrated-report  (calibration_report alone)
#   4. region_type = window                         -> OUTDIR/window           (as runs/example/snv-set-window)
#   5. region_type = coding                         -> OUTDIR/coding           (as runs/example/snv-set-coding)
#   6. region_type = coding, staar_enabled = true   -> OUTDIR/coding-staar     (mechanics only)
#   8. Rscript cloud/terra/check-wdl-example.R --wdl OUTDIR --reference REFERENCE [--tolerance T]
#
# Usage (from the GLOWanalyses root; miniwdl installed, Docker or Apptainer available):
#   bash cloud/terra/run-example-wdl.sh [OUTDIR]            # default OUTDIR = runs/example-wdl
# Environment:
#   MINIWDL                 the miniwdl executable (default: miniwdl on PATH)
#   GLOW_WDL_IMAGE          the image to run (default: the one glow_scan.wdl names)
#   GLOW_WDL_REFERENCE      a GLOWanalyses tree where run-example.sh has run (default: this tree)
#   MINIWDL__SCHEDULER__CONTAINER_BACKEND=singularity   for Apptainer instead of Docker (see
#                           cloud/README.md for the two run-option settings clusters need)
#   GLOW_WDL_SKIP_CHECK=1   stop after the seven runs (the comparison is then run apart, for
#                           example inside the image when the host has no R)
#   GLOW_WDL_TOLERANCE      a relative tolerance for the comparison's two full-precision
#                           artifacts (ld_scores.rds, the calibrated single-variant table), for a
#                           reference computed in another environment; default 0 (exact)
# Exit status: non-zero if a run or the comparison fails.
#
# Derived from the Terra workflows written by Alberto Brusati (Landers lab, UMass Chan Medical
# School), 2026.
set -euo pipefail
cd "$(dirname "$0")/../.."    # the GLOWanalyses root
A="$(pwd)"
OUT="${1:-runs/example-wdl}"
MINIWDL="${MINIWDL:-miniwdl}"
REF="${GLOW_WDL_REFERENCE:-$A}"
WDL="cloud/terra/glow_scan.wdl"
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

img_arg=()
if [[ -n "${GLOW_WDL_IMAGE:-}" ]]; then img_arg=("docker_image=${GLOW_WDL_IMAGE}"); fi

run_wdl() {   # run_wdl <subdir> <inputs.json> [extra input overrides...]
  local sub="$1"; local inputs="$2"; shift 2
  local dir="${OUT:?}/${sub:?}"          # both set, or the shell stops here
  if [[ -d "$dir" ]]; then rm -r -- "$dir"; fi
  mkdir -p "$dir"
  echo "=== miniwdl run: $sub ==="
  "$MINIWDL" run "$WDL" "$@" "${img_arg[@]}" -i "$inputs" -d "$dir/." --verbose > "$OUT/$sub.log" 2>&1 \
    || { echo "ERROR: miniwdl run $sub failed; see $OUT/$sub.log" >&2; tail -40 "$OUT/$sub.log" >&2; exit 1; }
  grep -E "workflow (start|done)|task .* done|error" "$OUT/$sub.log" | tail -12 || true
  test -f "$dir/outputs.json" || { echo "ERROR: $sub produced no outputs.json" >&2; exit 1; }
}

run_wdl single-variant cloud/terra/glow_scan.example.single_variant.inputs.json
run_wdl gene           cloud/terra/glow_scan.example.gene.inputs.json
SV_SCAN=$(find "$OUT/single-variant/out" -name marginal_all.csv | head -1)
SV_LD=$(find "$OUT/single-variant/out" -name ld_scores.rds | head -1)
test -n "$SV_SCAN" && test -n "$SV_LD" || { echo "ERROR: the single-variant run left no marginal_all.csv or ld_scores.rds" >&2; exit 1; }
run_wdl gene-calibrated cloud/terra/glow_scan.example.gene.inputs.json "single_variant_scan=$SV_SCAN" "ld_scores=$SV_LD"
SV_REPORT=$(find "$OUT/single-variant/out" -name single_variant_calibration.tsv | head -1)
test -n "$SV_REPORT" || { echo "ERROR: the single-variant run left no single_variant_calibration.tsv" >&2; exit 1; }
run_wdl gene-calibrated-report cloud/terra/glow_scan.example.gene.inputs.json "calibration_report=$SV_REPORT"
run_wdl window         cloud/terra/glow_scan.example.window.inputs.json
run_wdl coding         cloud/terra/glow_scan.example.coding.inputs.json
run_wdl coding-staar   cloud/terra/glow_scan.example.coding.inputs.json "staar_enabled=true"

if [[ "${GLOW_WDL_SKIP_CHECK:-0}" == "1" ]]; then
  echo "=== seven runs done under $OUT; the comparison is skipped (GLOW_WDL_SKIP_CHECK=1) ==="
  exit 0
fi
echo "=== compare with the run-example.sh outputs in $REF ==="
Rscript cloud/terra/check-wdl-example.R --wdl "$OUT" --reference "$REF" ${GLOW_WDL_TOLERANCE:+--tolerance "$GLOW_WDL_TOLERANCE"}
