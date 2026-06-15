#!/usr/bin/env bash
# ============================================================================
# data-example/build-intermediates.sh - build the TRACKED pre-built intermediates
#
# Runs the early GLOWanalyses templates on the synthetic example cohort to produce
# the pre-built, TRACKED derived assets so 03-snv-set / 04-summary (and the 01/02
# consumers) can run STANDALONE on a fresh clone:
#   - example_gds/chr{21,22}.gds                      (00/plink-to-gds)
#   - example_gds_favor/flexible/{gds,csv}/chr{21,22}_favor.{gds,csv}  (00/annotate-favor)
#       the aGDS is the scannable STAARpipeline sub-node-folder format that
#       annotate_favor now writes (carries the 16 PI features + the 3 coding nodes)
#   - example_pcs/pcs.rds                             (00/compute-pcs)
#   - example_pheno/data-prep_pheno_covar.rds         (00/assemble-pheno-covar)
#   - cohort/b_model.rds                              (01-training/b -> copied up)
#   - cohort/pi_models/                               (01-training/pi -> copied up)
#
# Run AFTER make-example-cohort.R. Deterministic + idempotent. With GLOWr +
# GLOWpipeline installed (the script cd's to the GLOWanalyses root, so it runs from
# anywhere):
#   conda activate r_env
#   bash data-example/build-intermediates.sh
set -euo pipefail

# Anchor to the GLOWanalyses root (this script lives in data-example/).
cd "$(dirname "$0")/.."
ROOT="data-example/cohort"
DP_CFG="runs/example/data-prep/config.R"
TR_CFG="runs/example/training/config.R"
A="."

if [[ ! -f "$ROOT/pheno.csv" ]]; then
  echo "ERROR: raw cohort not found. Run make-example-cohort.R first." >&2
  exit 1
fi

echo "=================================================================="
echo "  Building example pre-built intermediates"
echo "=================================================================="

echo "--- [1/6] 00/plink-to-gds ---"
Rscript "$A/00-data-prep/plink-to-gds.R" --config "$DP_CFG"

echo "--- [2/6] 00/annotate-favor (produces the scannable aGDS + its CSV) ---"
# annotate_favor writes /annotation/info/FunctionalAnnotation as a STAARpipeline
# sub-node FOLDER (one native-typed sub-node per feature), so the produced aGDS is
# directly scannable by 02/03 AND carries the 3 string coding nodes the coding
# masks read. The example data-prep config sets favor_features = 16 PI + 3 coding.
Rscript "$A/00-data-prep/annotate-favor.R" --config "$DP_CFG"

echo "--- [3/6] 00/compute-pcs ---"
Rscript "$A/00-data-prep/compute-pcs.R" --config "$DP_CFG"

echo "--- [4/6] 00/assemble-pheno-covar ---"
Rscript "$A/00-data-prep/assemble-pheno-covar.R" --config "$DP_CFG"

echo "--- [5/6] 01-training/b (prepare -> estimate) ---"
Rscript "$A/01-training/b/prepare-b-training-data.R" --config "$TR_CFG"
Rscript "$A/01-training/b/estimate-b.R"              --config "$TR_CFG"

echo "--- [6/6] 01-training/pi (prepare -> annotate -> train -> evaluate) ---"
Rscript "$A/01-training/pi/prepare-pi-case.R"  --config "$TR_CFG"
Rscript "$A/01-training/pi/annotate-pi-case.R" --config "$TR_CFG"
Rscript "$A/01-training/pi/train-pi.R"         --config "$TR_CFG"
Rscript "$A/01-training/pi/evaluate-pi.R"      --config "$TR_CFG"

echo "--- Copy the pre-built derived assets up to the cohort root (TRACKED) ---"
# The aGDS already lives at its tracked path (00/annotate-favor output). Copy the B
# model, the PI ensemble, and the pheno bundle to the cohort root so the ANALYSIS
# stages (02-single-variant, 03-snv-set, 04-summary) run STANDALONE on a fresh clone,
# with none of the git-ignored 00/01 lineage dirs present.
TR_OUT="$A/runs/example/training/outputs"
cp -f "$TR_OUT/b_model.rds" "$ROOT/b_model.rds"
rm -rf "$ROOT/pi_models"
cp -r "$TR_OUT/pi_models" "$ROOT/pi_models"
cp -f "$ROOT/example_pheno/data-prep_pheno_covar.rds" "$ROOT/pheno_covar.rds"
echo "  cohort/b_model.rds + cohort/pi_models/ + cohort/pheno_covar.rds updated."

echo "=================================================================="
echo "  Intermediates built. Tracked pre-built assets:"
echo "    $ROOT/example_gds_favor/flexible/gds/chr{21,22}_favor.gds"
echo "    $ROOT/b_model.rds ; $ROOT/pi_models/"
echo "=================================================================="
