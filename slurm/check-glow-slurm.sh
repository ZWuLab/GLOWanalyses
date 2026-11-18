#!/usr/bin/env bash
# ==============================================================================
# File: slurm/check-glow-slurm.sh
# Purpose: Check the SLURM mode of the `glow` command on a cluster, with the shipped
#   example. It submits real jobs (about twenty, each of a few seconds) and checks:
#
#   (a) one call finishes the example's gene run;
#   (b) a second call submits nothing and reports every unit as finished;
#   (c) after a deleted unit result, a call submits that unit and the parts after
#       it only;
#   (d) the window run finishes from one call with --wait;
#   (e) a unit that fails is submitted again by the next call, together with the
#       units that read it;
#   (f) the tables of the gene and window runs equal those of a local run of the
#       same steps, within a relative tolerance of 1e-8, because the jobs may run on
#       other processors than the local run;
#   (g) a job whose configuration was changed after its submission does not run,
#       and says so in its log;
#   (h) a call with --wait stops with an error when a job it follows fails, and the
#       same command then finishes the run.
#
# Usage (from a shell where `sbatch` works and the R environment is active):
#   bash slurm/check-glow-slurm.sh <work folder on storage that the compute nodes see> [sbatch options]
#
# The work folder has to be empty or absent. The second argument is passed to
# every sbatch call through run.slurm.options, for example "--partition=short".
# The jobs set up R by activating the conda environment named in GLOW_CONDA_ENV
# (default r_env). The check takes five to ten minutes on an idle cluster. Each
# check prints PASS or FAIL, and the last line is "Verdict for the SLURM mode: PASS"
# or "... FAIL". The exit status is 1 on a failure. Logs are kept in the work
# folder.
set -uo pipefail
cd "$(dirname "$0")/.."
GA="$PWD"
W="${1:?give a work folder on storage that the compute nodes see}"
OPTIONS="${2:-}"
mkdir -p "$W"
if [[ -n "$(ls -A "$W")" ]]; then echo "ERROR: the work folder $W is not empty." >&2; exit 1; fi
W="$(cd "$W" && pwd)"
command -v sbatch >/dev/null 2>&1 || { echo "ERROR: sbatch is not on PATH." >&2; exit 1; }
echo "Work folder: $W"

# ---- the configurations --------------------------------------------------------
slurm_block() { cat <<EOF2
run:
  mode: slurm
  slurm:
    options: "$OPTIONS"
    conda_env: "${GLOW_CONDA_ENV:-r_env}"
    poll_seconds: 15
    resources:
      default: {mem: 4G, time: "00:20:00"}
EOF2
}
{ echo "extends: $GA/example.yml"; echo "output_root: $W/runs"; slurm_block; echo "snv_set:"; echo "  region_type: [gene]"; } > "$W/gene.yml"
{ echo "extends: $GA/example.yml"; echo "output_root: $W/runs"; slurm_block; echo "snv_set:"; echo "  region_type: [window]"; } > "$W/window.yml"
# a run of its own analysis in which the job of chr22 fails once: the setup line ends the
# first job of that unit with an error and leaves a flag file, so the next one runs
# (it extends gene.yml, so it reads the pre-built files of the analysis baseline, and it sets
# top_k, so its gene run is its own and not a reuse of the one above)
{ echo "extends: gene.yml"; echo "analysis: fails_once"; slurm_block
  echo "    setup:"
  echo "      - 'if [[ \"\$SLURM_JOB_NAME\" == *scan-gene-chr22-* && ! -e $W/failed-once ]]; then touch $W/failed-once; echo \"planted failure\"; exit 7; fi'"
  echo "snv_set:"; echo "  top_k: 10"; } > "$W/fails_once.yml"
{ echo "extends: $GA/example.yml"; echo "output_root: $W/local"; echo "snv_set:"; echo "  region_type: [gene, window]"; } > "$W/local.yml"
# the same planted failure in an analysis of its own, for a call with --wait
{ echo "extends: gene.yml"; echo "analysis: fails_wait"; slurm_block
  echo "    setup:"
  echo "      - 'if [[ \"\$SLURM_JOB_NAME\" == *scan-gene-chr22-* && ! -e $W/failed-wait ]]; then touch $W/failed-wait; echo \"planted failure\"; exit 7; fi'"
  echo "snv_set:"; echo "  top_k: 12"; } > "$W/fails_wait.yml"
# a run of its own analysis whose jobs start a minute after their submission, so that the
# configuration can be changed in between (--begin is an option of the site's kind, which
# run.slurm.options passes on)
{ echo "extends: gene.yml"; echo "analysis: changed"; OPTIONS="$OPTIONS --begin=now+75" slurm_block
  echo "snv_set:"; echo "  top_k: 11"; } > "$W/changed.yml"

# ---- helpers -------------------------------------------------------------------
N=0; FAILED=0; LAST_LOG=""
check() {   # check <name> <ok|stop> <pattern or -> <command...>
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
  [[ "$verdict" == PASS ]] || { echo "      log: $log"; tail -n 6 "$log" | sed 's/^/      | /'; }
}
absent() { ! grep -qE -- "$1" "$2"; }
# wait until none of this check's jobs is in the queue (their names carry the analysis)
drain() { local tag="$1" i; for i in $(seq 1 120); do
            [[ "$(squeue -h -u "$USER" -o '%j' | grep -c -- "^glow-example-$tag-")" == "0" ]] && return 0; sleep 10; done
          echo "the jobs of $tag did not leave the queue in 20 minutes"; return 1; }

echo "--- (a) one call finishes the gene run ---"
check "a: the call submits the five units" ok "plots gene +submitted as job [0-9]+, afterok:[0-9]+" ./glow run snv_set --config "$W/gene.yml"
check "a: the scan units wait for the preparation" ok "scan gene chr21 +submitted as job [0-9]+, afterok:[0-9]+" cat "$LAST_LOG"
check "a: the jobs leave the queue" ok - drain baseline
echo "--- (b) a second call submits nothing ---"
check "b: every unit is finished" ok "every unit is finished" ./glow run snv_set --config "$W/gene.yml"
check "b: nothing is submitted" ok - absent "submitted as job" "$LAST_LOG"

echo "--- (c) a deleted result: that unit and the parts after it only ---"
rm -f "$W/runs/example/baseline/snv_set/gene/results/glow_chr22.csv"
check "c: the call submits chr22, the aggregation and the plots" ok "scan gene chr22 +submitted as job" ./glow run snv_set --config "$W/gene.yml"
C_LOG="$LAST_LOG"
check "c: the aggregation follows the scan" ok "aggregate gene +submitted as job [0-9]+, afterok:[0-9]+" cat "$C_LOG"
check "c: the plots follow the aggregation" ok "plots gene +submitted as job [0-9]+, afterok:[0-9]+" cat "$C_LOG"
check "c: the preparation is left alone" ok "prepare gene +complete" cat "$C_LOG"
check "c: the other chromosome is left alone" ok "scan gene chr21 +complete" cat "$C_LOG"
check "c: the jobs leave the queue" ok - drain baseline
check "c: the result file is back" ok - test -s "$W/runs/example/baseline/snv_set/gene/results/glow_chr22.csv"
check "c: every unit is finished again" ok "every unit is finished" ./glow run snv_set --config "$W/gene.yml"

echo "--- (d) the window run from one call with --wait ---"
check "d: the call stays until the run is finished" ok "every unit is finished" ./glow run snv_set --config "$W/window.yml" --wait
check "d: the chunks were submitted by the same call" ok "scan window chunk0001 +submitted as job" cat "$LAST_LOG"
check "d: the summary waited for the chunks" ok "aggregate window +waits for the units of scan" cat "$LAST_LOG"

echo "--- (e) a failed unit is submitted again with what reads it ---"
check "e: the first call submits the run" ok "scan gene chr22 +submitted as job" ./glow run snv_set --config "$W/fails_once.yml"
check "e: the jobs leave the queue" ok - drain fails_once
check "e: the unit failed once" ok - test -e "$W/failed-once"
check "e: the next call submits chr22 again, and the summary" ok "scan gene chr22 +submitted as job" ./glow run snv_set --config "$W/fails_once.yml"
E_LOG="$LAST_LOG"
check "e: the aggregation is submitted again" ok "aggregate gene +submitted as job [0-9]+, afterok:[0-9]+" cat "$E_LOG"
check "e: what had finished is left alone" ok "scan gene chr21 +complete" cat "$E_LOG"
check "e: the jobs leave the queue" ok - drain fails_once
check "e: every unit is finished" ok "every unit is finished" ./glow run snv_set --config "$W/fails_once.yml"

echo "--- (h) --wait stops when a job it follows fails, and the same command then finishes the run ---"
check "h: the call stops and names the unit" stop "left the queue without finishing: .*scan gene chr22" ./glow run snv_set --config "$W/fails_wait.yml" --wait
H_LOG="$LAST_LOG"
check "h: it submitted each of the five units once" ok "^5$" grep -c "submitted as job" "$H_LOG"
check "h: the same command finishes the run" ok "every unit is finished" ./glow run snv_set --config "$W/fails_wait.yml" --wait
check "h: what had finished was left alone" ok "scan gene chr21 +complete" cat "$LAST_LOG"

echo "--- (g) a job does not run with a configuration that changed after its submission ---"
check "g: the call submits the run, to start in a minute" ok "prepare gene +submitted as job" ./glow run snv_set --config "$W/changed.yml"
echo "# a line added after the submission" >> "$W/changed.yml"
check "g: the jobs leave the queue" ok - drain changed
check "g: the job says why it did not run" ok "a file of the configuration changed after this job was submitted" \
  bash -c "cat $W/runs/example/changed/snv_set/gene/logs/slurm/*prepare*.log"
check "g: no unit of the run has a record" stop - bash -c "ls $W/runs/example/changed/snv_set/gene/records/*.yml"
check "g: the next call would submit the run again" ok "prepare gene +would submit" ./glow run snv_set --config "$W/changed.yml" --dry-run

echo "--- (f) the tables equal those of a local run ---"
check "f: the local run of the two region types" ok - ./glow run snv_set --config "$W/local.yml"
check "f: the comparison within 1e-8" ok "Verdict for the comparison: PASS" Rscript data-example/compare-example-runs.R \
  --a "$GA" --a-form yaml --a-root "$W/local" --b "$GA" --b-form yaml --b-root "$W/runs" --groups gene,window --tolerance 1e-8 --out "$W/comparison.md"

echo
echo "Checks: $N, failed: $FAILED. Logs: $W"
if [[ $FAILED -eq 0 ]]; then echo "Verdict for the SLURM mode: PASS"; else echo "Verdict for the SLURM mode: FAIL"; exit 1; fi
