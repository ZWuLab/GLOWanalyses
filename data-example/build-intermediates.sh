#!/usr/bin/env bash
# ============================================================================
# data-example/build-intermediates.sh - build the pre-built files of the example
#
# The example ships four kinds of files that a full run would produce, so that a
# new user can start at the analysis steps (data-example/seed-manifest.yml lists
# them): the annotated GDS files, the phenotype bundle, the B model and the PI
# models. This script builds them again from the raw synthetic cohort and copies
# them into data-example/cohort/, where they are tracked:
#
#   ./glow run data_prep --force      PLINK to GDS, FAVOR annotation, PCs, phenotype bundle
#   ./glow run training --force       the B model and the PI models
#   cohort/example_gds_favor/flexible/gds/chr{21,22}_favor.gds   <- runs/example/data/gds_favor/flexible/gds/
#   cohort/pheno_covar.rds                                       <- runs/example/baseline/data_prep/
#   cohort/b_model.rds                                           <- runs/example/baseline/training/b/
#   cohort/pi_models/ and cohort/pi_models.tar.gz                <- runs/example/baseline/training/pi/pi_models/
#
# Run it after make-example-cohort.R, when the generator or a setting of the two
# steps changed. It is a maintainer's tool: a user never needs it. With GLOWr and
# GLOWpipeline installed (the script runs from any directory):
#   conda activate r_env
#   bash data-example/build-intermediates.sh
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="data-example/cohort"

if [[ ! -f "$ROOT/pheno.csv" ]]; then
  echo "ERROR: raw cohort not found. Run make-example-cohort.R first." >&2
  exit 1
fi

echo "=================================================================="
echo "  Building the pre-built files of the example"
echo "=================================================================="
# --wait: with run.mode slurm the command stays until the jobs have finished, and stops when one failed.
./glow run data_prep --force --wait
./glow run training --force --wait

echo "--- copy the rebuilt files into $ROOT (tracked), and fix the settings they were built with ---"
# The script takes the folders from the layout of example.yml, copies a file only when a run of
# example.yml as it is now finished its unit, and then writes data-example/seed-settings.yml:
# the seed records of a clone state these settings, whatever example.yml says later.
Rscript data-example/write-seed-settings.R --collect
# the PI models as the one tarball that the Terra workflow takes
tar -czf "$ROOT/pi_models.tar.gz" -C "$ROOT" pi_models

echo "=================================================================="
echo "  Done. Tracked pre-built files:"
echo "    $ROOT/example_gds_favor/flexible/gds/chr{21,22}_favor.gds"
echo "    $ROOT/pheno_covar.rds ; $ROOT/b_model.rds ; $ROOT/pi_models/ ; $ROOT/pi_models.tar.gz"
echo "    data-example/seed-settings.yml (the settings they were built with)"
echo "  The seed records of a clone are written from these files on its first run."
echo "=================================================================="
