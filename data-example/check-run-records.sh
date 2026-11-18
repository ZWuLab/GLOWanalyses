#!/usr/bin/env bash
# ==============================================================================
# File: data-example/check-run-records.sh
# Purpose: Check the records of a run on the shipped example: a finished unit is
#   skipped, a derived analysis reuses what did not change and is told what to run
#   again, a run that would overwrite results made with other settings is refused,
#   a forced rerun makes what read its output incomplete, a deleted result makes
#   its unit and only its unit run again, an input from the example data stops a
#   run of another dataset, a finished folder of the parent is reused whole and an
#   unfinished one is not, and a setting on the command line is refused.
#
#   The script runs the stage scripts on configurations that extend example.yml and
#   write to a work folder of their own, so runs/example/ is not touched.
#
# Usage (from the GLOWanalyses directory, or anywhere):
#   bash data-example/check-run-records.sh [work folder]
#
# The work folder defaults to a new temporary folder. It has to be empty or absent:
# the checks start from a dataset without records. Each check prints PASS or FAIL,
# and the last line is "Verdict for the run records: PASS" or "... FAIL". The exit
# status is 1 on a failure. The output of every command is kept in the work folder
# (NN-<check>.log).
set -uo pipefail
cd "$(dirname "$0")/.."
GA="$PWD"
W="${1:-$(mktemp -d)}"
mkdir -p "$W"
if [[ -n "$(ls -A "$W")" ]]; then echo "ERROR: the work folder $W is not empty." >&2; exit 1; fi
W="$(cd "$W" && pwd)"
echo "Work folder: $W"

# ---- the configurations of the checks ------------------------------------------
cat > "$W/base.yml" <<EOF
# the shipped example, one region type, with its outputs in the work folder
extends: $GA/example.yml
output_root: $W/runs
snv_set:
  region_type: [gene]
EOF
cat > "$W/change_pi.yml" <<EOF
# a derived analysis: three of the sixteen PI features
extends: base.yml
analysis: change_pi
shared:
  pi_features: [cadd_phred, linsight, fathmm_xf]
EOF
cat > "$W/changed_ld.yml" <<EOF
# the same analysis as base.yml, with one changed setting
extends: base.yml
snv_set:
  ld_threshold: 0.8
EOF
cat > "$W/same.yml" <<EOF
# a derived analysis that changes nothing the gene run reads
extends: base.yml
analysis: same
EOF
cat > "$W/parent2.yml" <<EOF
# a second analysis, of which only the preparation of the gene run is finished
extends: base.yml
analysis: parent2
snv_set:
  ld_threshold: 0.7
EOF
cat > "$W/child2.yml" <<EOF
# a derived analysis of the unfinished one, with the same settings
extends: parent2.yml
analysis: child2
EOF
cat > "$W/not_example.yml" <<EOF
# another dataset whose raw inputs were left at their defaults, the example's files
dataset: mine
output_root: $W/runs
EOF

# ---- how a check is run and judged ---------------------------------------------
N=0; FAILED=0; LAST_LOG=""
# check <name> <ok|stop> <pattern or -> <command...>
#   ok    the command has to succeed, and its output has to match the pattern
#   stop  the command has to fail, and its output has to match the pattern
#   A pattern of "-" asks for nothing. LAST_LOG is the log of the check before, so a
#   second check can look at the same output.
check() {
  local name="$1" want="$2" pattern="$3"; shift 3
  N=$((N + 1)); local log; log="$W/$(printf '%02d' "$N")-$(echo "$name" | tr -c 'A-Za-z0-9\n' '-').log"
  "$@" > "$log" 2>&1; local status=$?
  LAST_LOG="$log"
  local verdict=PASS
  if [[ "$want" == ok && $status -ne 0 ]]; then verdict="FAIL (exit $status, expected success)"; fi
  if [[ "$want" == stop && $status -eq 0 ]]; then verdict="FAIL (succeeded, expected a stop)"; fi
  if [[ "$pattern" != "-" ]] && ! grep -qE -- "$pattern" "$log"; then verdict="FAIL (no match for: $pattern)"; fi
  [[ "$verdict" == PASS ]] || FAILED=$((FAILED + 1))
  printf '  [%s] %s\n' "$verdict" "$name"
  [[ "$verdict" == PASS ]] || { echo "      log: $log"; tail -n 5 "$log" | sed 's/^/      | /'; }
}
gene() { Rscript 03-snv-set/run-gene.R --config "$1" --chr "$2" "${@:3}"; }
# the checksum of every result table of an analysis's gene run
results_sum() { (cd "$W/runs/example/$1/snv_set/gene" && find results aggregated shared -type f | sort | xargs md5sum | md5sum | cut -d' ' -f1); }

echo "--- a first run of the gene scan on the example ---"
check "prepare" ok - Rscript 03-snv-set/prepare.R --config "$W/base.yml"
check "scan chr21" ok - gene "$W/base.yml" 21
check "scan chr22" ok - gene "$W/base.yml" 22
check "aggregate" ok - Rscript 04-summary/aggregate.R --config "$W/base.yml"
check "plots" ok - Rscript 04-summary/plots.R --config "$W/base.yml"
check "a finished unit is skipped" ok "scan chr21: complete .*Nothing to do" gene "$W/base.yml" 21

echo "--- (g) a derived analysis with the same settings reuses the finished folder of its parent, whole ---"
check "g: the preparation is reused" ok "prepare: reused from the analysis .baseline." Rscript 03-snv-set/prepare.R --config "$W/same.yml"
check "g: a scan unit is reused" ok "scan chr22: reused from the analysis .baseline." gene "$W/same.yml" 22
check "g: the aggregation is reused" ok "aggregate: reused from the analysis .baseline." Rscript 04-summary/aggregate.R --config "$W/same.yml"
check "g: nothing is written in the derived analysis" stop - test -e "$W/runs/example/same/snv_set"

echo "--- (h) a parent with a finished preparation and no scan is not a finished folder ---"
check "h: the parent's preparation runs" ok - Rscript 03-snv-set/prepare.R --config "$W/parent2.yml"
check "h: the scan of the derived analysis is not reused from results that do not exist" stop "snv_set, part prepare" gene "$W/child2.yml" 21
check "h: its message does not say reused" stop - grep -q "reused from" "$LAST_LOG"
check "h: the derived analysis runs its own preparation" ok "Stage 1 complete" Rscript 03-snv-set/prepare.R --config "$W/child2.yml"
check "h: and then its own scan" ok "Done gene chr21" gene "$W/child2.yml" 21

echo "--- (i) a setting on the command line of a YAML run is refused ---"
check "i: --max-genes stops the scan and names the setting" stop "--max-genes is not an option of a YAML configuration. Set .snv_set.max_genes." \
  gene "$W/child2.yml" 22 --max-genes 1
check "i: no record was written for the refused unit" stop - test -e "$W/runs/example/child2/snv_set/gene/records/scan.chr22.yml"
check "i: and no lock is held" stop - bash -c "compgen -G \"$W/runs/example/child2/snv_set/gene/records/scan.chr22.lock.*/owner\" > /dev/null"

echo "--- (j) a deleted file inside a recorded folder makes its unit run again ---"
rm -f "$W/runs/example/parent2/snv_set/gene/shared/01-config_used.rds"
check "j: the preparation is not complete without its settings record" ok "Stage 1 complete" Rscript 03-snv-set/prepare.R --config "$W/parent2.yml"
check "j: the file is back" ok - test -s "$W/runs/example/parent2/snv_set/gene/shared/01-config_used.rds"

echo "--- (a) a derived analysis with other PI features stops at the SNV-set step and names the training ---"
check "a: the preparation stops and names training" stop "training, part pi_train .*glow run training" \
  Rscript 03-snv-set/prepare.R --config "$W/change_pi.yml"
check "a: the message names the setting" ok "setting shared.pi_features" cat "$LAST_LOG"

echo "--- (b) after the PI training it runs, on the supplied B model and phenotype bundle of the parent ---"
check "b: the PI cases are prepared" ok - Rscript 01-training/pi/prepare-pi-case.R --config "$W/change_pi.yml"
check "b: the PI cases are annotated" ok - Rscript 01-training/pi/annotate-pi-case.R --config "$W/change_pi.yml"
check "b: new PI models are trained" ok - Rscript 01-training/pi/train-pi.R --config "$W/change_pi.yml"
check "b: the new models are in the derived analysis" ok "model_1.rds" ls "$W/runs/example/change_pi/training/pi/pi_models"
check "b: no B model is written in the derived analysis" stop - test -e "$W/runs/example/change_pi/training/b/b_model.rds"
check "b: the preparation now runs" ok - Rscript 03-snv-set/prepare.R --config "$W/change_pi.yml"
RESOLVED="$W/runs/example/change_pi/snv_set/gene/config-resolved.yml"
check "b: the record states the reuse of the B model" ok "training/b_estimate: reused from baseline" cat "$RESOLVED"
check "b: the record states the reuse of the bundle" ok "data_prep/assemble_pheno: reused from baseline" cat "$RESOLVED"
check "b: the B model path is the parent's" ok "b_model_path: .*/runs/example/baseline/training/b/b_model.rds" cat "$RESOLVED"
check "b: the PI models are the analysis's own" ok "pi_model_dir: .*/runs/example/change_pi/training/pi/pi_models" cat "$RESOLVED"
check "b: the scan of the derived analysis runs" ok - gene "$W/change_pi.yml" 21
check "b: the record still states the reuse after the scan wrote it again" ok "training/b_estimate: reused from baseline" cat "$RESOLVED"

echo "--- (c) a changed setting aimed at a finished folder is refused without --force ---"
SUM_BEFORE="$(results_sum baseline)"
check "c: the preparation is refused" stop "setting snv_set.ld_threshold: recorded 0.9.*now 0.8" \
  Rscript 03-snv-set/prepare.R --config "$W/changed_ld.yml"
check "c: the refusal gives the two ways forward" ok "a new .analysis. name.*--force" cat "$LAST_LOG"
check "c: the finished unit is still complete" ok "scan chr22: complete" gene "$W/base.yml" 22
check "c: every result file has the checksum it had before the refusal" ok - test "$SUM_BEFORE" = "$(results_sum baseline)"

echo "--- (d) a change of ld_threshold makes the scan units of a finished run incomplete ---"
check "d: the forced preparation runs" ok - Rscript 03-snv-set/prepare.R --config "$W/changed_ld.yml" --force
check "d: a scan unit of the old run is not complete any more" stop "made with other dependencies" gene "$W/changed_ld.yml" 21
check "d: under the old configuration the preparation is what is missing" stop "snv_set, part prepare" gene "$W/base.yml" 22
check "d: the forced scan of chr21 runs" ok - gene "$W/changed_ld.yml" 21 --force
check "d: the forced scan of chr22 runs" ok - gene "$W/changed_ld.yml" 22 --force

echo "--- (e) a deleted result file makes its unit, and only its unit, run again ---"
rm -f "$W/runs/example/baseline/snv_set/gene/results/glow_chr22.csv"
check "e: the other unit is still complete" ok "scan chr21: complete .*Nothing to do" gene "$W/changed_ld.yml" 21
check "e: the unit of the deleted file runs again" ok "Done gene chr22" gene "$W/changed_ld.yml" 22
check "e: the result file is back" ok - test -s "$W/runs/example/baseline/snv_set/gene/results/glow_chr22.csv"

echo "--- (f) an input from the example data stops a run of another dataset ---"
check "f: the run stops and names the setting" stop "belong to the shipped example" \
  Rscript 01-training/b/prepare-b-training-data.R --config "$W/not_example.yml"
check "f: the setting is training.b_raw_path" ok "training.b_raw_path" cat "$LAST_LOG"

echo
echo "Checks: $N, failed: $FAILED. Logs: $W"
if [[ $FAILED -eq 0 ]]; then echo "Verdict for the run records: PASS"; else echo "Verdict for the run records: FAIL"; exit 1; fi
