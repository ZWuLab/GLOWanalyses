#!/usr/bin/env bash
# ============================================================================
# File: data-example/check-seed-collector.sh
# Purpose: The regression check of the collector of the example's pre-built files
#   (data-example/write-seed-settings.R --collect). The collector may put files into
#   the tracked cohort only when they are the result of finished units of example.yml
#   as it is, and the files and data-example/seed-settings.yml have to be one set.
#   This script disturbs the collector in the ways a second process could, one at a
#   time, and checks that nothing is collected and that the tracked files keep their
#   bytes. It then interrupts a collection and checks that neither --check nor a run
#   of the example accepts what is left.
#
# Usage:
#   bash data-example/check-seed-collector.sh <work folder>
#
#   The script copies this GLOWanalyses tree into <work folder>/ga and works there. It
#   never changes the tree it is called from. <work folder> must not exist or must be
#   empty. The copy runs data-example/build-intermediates.sh first (the data
#   preparation and the training of the example, a few minutes), with the R that is
#   on PATH.
#
# Output: one line per check, `Checks: <n>, failed: <m>.` and
#   `Verdict for the seed collector: PASS|FAIL`. Exit status 1 on FAIL. The logs of
#   every step are in <work folder>.
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
W="${1:-}"
if [[ -z "$W" ]]; then echo "Usage: bash data-example/check-seed-collector.sh <work folder>" >&2; exit 2; fi
mkdir -p "$W"; W="$(cd "$W" && pwd)"
if [[ -n "$(ls -A "$W")" ]]; then echo "ERROR: the work folder $W is not empty." >&2; exit 2; fi
case "$W/" in "$HERE"/*) echo "ERROR: the work folder must lie outside the GLOWanalyses tree." >&2; exit 2;; esac

# the copy: the tracked files, without runs and job folders
mkdir "$W/ga"
(cd "$HERE" && tar -cf - --exclude='./runs' --exclude='./.git' --exclude='outputs' --exclude='jobs' .) | (cd "$W/ga" && tar -xf -)
cd "$W/ga"

bash data-example/build-intermediates.sh > "$W/builder.log" 2>&1; st=$?
echo "builder: exit $st ($(grep -c '^collected ' "$W/builder.log") file(s) or folder(s) collected)"
if [[ $st -ne 0 ]]; then echo "ERROR: the builder failed; see $W/builder.log" >&2; exit 1; fi

MARKER=data-example/seed-collection.incomplete
# the tracked pre-built files and the settings file, as one checksum
sums() { find data-example/cohort/pheno_covar.rds data-example/cohort/b_model.rds data-example/cohort/pi_models data-example/cohort/example_gds_favor data-example/seed-settings.yml -type f | sort | xargs md5sum | md5sum | cut -d' ' -f1; }
BEFORE=$(sums); N=0; BAD=0
judge() { N=$((N + 1)); if [[ "$2" == yes ]]; then echo "  [PASS] $1"; else echo "  [FAIL] $1"; BAD=$((BAD + 1)); fi; }
same() { [[ "$(sums)" == "$BEFORE" && ! -e $MARKER ]] && echo yes || echo no; }
probe() { Rscript data-example/seed-collector-probe.R "$1" > "$W/probe-$1.log" 2>&1; }

# 1. example.yml is edited while the files are copied
cp example.yml "$W/example.yml.kept"
probe edit; st=$?
mv -f "$W/example.yml.kept" example.yml
judge "an edit of example.yml during the copy stops the collector" "$([[ $st -ne 0 ]] && grep -q 'A file of the configuration changed while the files were copied: example.yml' "$W/probe-edit.log" && echo yes || echo no)"
judge "and the tracked files and the settings file are unchanged, with no marker" "$(same)"

# 2. a collected unit is produced again while the files are copied
probe rerun; st=$?
mv -f runs/example/baseline/training/b/records/b_estimate.yml.kept runs/example/baseline/training/b/records/b_estimate.yml 2>/dev/null
mv -f b_model.kept runs/example/baseline/training/b/b_model.rds 2>/dev/null
judge "a unit that loses its record during the copy stops the collector" "$([[ $st -ne 0 ]] && grep -q 'The unit training.b_estimate was run again, lost its completion record, or is no longer finished' "$W/probe-rerun.log" && echo yes || echo no)"
judge "and the tracked files and the settings file are unchanged, with no marker" "$(same)"
judge "and the bytes of the unfinished unit are not in the tracked cohort" "$(grep -q 'INCOMPLETE PRODUCER OUTPUT' data-example/cohort/b_model.rds && echo no || echo yes)"

# 3. a raw input is changed while the files are copied
cp data-example/cohort/pheno.csv "$W/pheno.csv.kept"
probe raw; st=$?
mv -f "$W/pheno.csv.kept" data-example/cohort/pheno.csv
judge "a raw input that changes during the copy stops the collector" "$([[ $st -ne 0 ]] && grep -Eq 'is no longer finished for its inputs|An input of the pre-built files changed while they were copied' "$W/probe-raw.log" && echo yes || echo no)"
judge "and the tracked files and the settings file are unchanged, with no marker" "$(same)"

# 4. a part before the collected files, itself not collected, is produced again during the copy
PCS_RECORD=$(find runs/example -path '*/records/compute_pcs*.yml' | head -1)
if [[ -n "$PCS_RECORD" ]]; then cp "$PCS_RECORD" "$W/compute_pcs.yml.kept"; fi
probe upstream; st=$?
if [[ -n "$PCS_RECORD" ]]; then mv -f "$W/compute_pcs.yml.kept" "$PCS_RECORD"; fi
judge "a part before the files that is produced again during the copy stops the collector" "$([[ $st -ne 0 ]] && grep -q 'is no longer finished for its inputs' "$W/probe-upstream.log" && echo yes || echo no)"
judge "and the tracked files and the settings file are unchanged, with no marker" "$(same)"

# 5. the two modes together, and an unknown argument
Rscript data-example/write-seed-settings.R --collect --check > "$W/both.log" 2>&1; st=$?
judge "--collect with --check stops and changes nothing" "$([[ $st -ne 0 && "$(same)" == yes ]] && grep -q 'not both' "$W/both.log" && echo yes || echo no)"
Rscript data-example/write-seed-settings.R --colect > "$W/unknown.log" 2>&1; st=$?
judge "an unknown argument stops and changes nothing" "$([[ $st -ne 0 && "$(same)" == yes ]] && grep -q 'Unknown argument' "$W/unknown.log" && echo yes || echo no)"

# 6. a collection that is interrupted after its first file
probe interrupt; st=$?
judge "an interrupted collection leaves its marker" "$([[ $st -ne 0 && -e $MARKER ]] && echo yes || echo no)"
Rscript data-example/write-seed-settings.R --check > "$W/check-interrupted.log" 2>&1; st=$?
judge "--check fails while the marker is there" "$([[ $st -ne 0 ]] && grep -q 'Verdict for the seed settings: FAIL' "$W/check-interrupted.log" && echo yes || echo no)"
Rscript data-example/write-seed-settings.R > "$W/plain-interrupted.log" 2>&1; st=$?
judge "the settings file is not written again over an interrupted collection" "$([[ $st -ne 0 ]] && grep -q 'began and did not end' "$W/plain-interrupted.log" && echo yes || echo no)"
./glow run single_variant --dry-run > "$W/run-interrupted.log" 2>&1; st=$?
judge "a run of the example does not use the files of an interrupted collection" "$([[ $st -ne 0 ]] && grep -q 'began and did not end' "$W/run-interrupted.log" && echo yes || echo no)"

# 7. an undisturbed collection repairs it and gives the same bytes
Rscript data-example/write-seed-settings.R --collect > "$W/collect-clean.log" 2>&1; st=$?
judge "an undisturbed collection succeeds, removes the marker and gives the same bytes" "$([[ $st -eq 0 && "$(same)" == yes ]] && echo yes || echo no)"
Rscript data-example/write-seed-settings.R --check > "$W/check-clean.log" 2>&1; st=$?
judge "--check passes afterward" "$([[ $st -eq 0 ]] && grep -q 'Verdict for the seed settings: PASS' "$W/check-clean.log" && echo yes || echo no)"

# 8. a pre-built file of another collection beside the settings file (no marker)
cp data-example/cohort/b_model.rds "$W/b_model.rds.kept"
printf 'another model' > data-example/cohort/b_model.rds
Rscript data-example/write-seed-settings.R --check > "$W/check-mixed.log" 2>&1; st=$?
judge "--check fails on a pre-built file that the settings file does not state" "$([[ $st -ne 0 ]] && grep -q 'Verdict for the seed settings: FAIL' "$W/check-mixed.log" && echo yes || echo no)"
./glow run single_variant --dry-run > "$W/run-mixed.log" 2>&1; st=$?
judge "and a run of the example does not use that set" "$([[ $st -ne 0 ]] && grep -q 'not the file that was collected: b_model.rds' "$W/run-mixed.log" && echo yes || echo no)"
mv -f "$W/b_model.rds.kept" data-example/cohort/b_model.rds
./glow run single_variant --dry-run > "$W/run-restored.log" 2>&1; st=$?
judge "with the file back, the run of the example reads its plan again" "$([[ $st -eq 0 && "$(same)" == yes ]] && echo yes || echo no)"

echo "Checks: $N, failed: $BAD."
if [[ $BAD -eq 0 ]]; then echo "Verdict for the seed collector: PASS"; else echo "Verdict for the seed collector: FAIL"; exit 1; fi
