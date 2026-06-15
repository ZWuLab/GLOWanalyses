# Synthetic example cohort (`data-example/`)

A small, **deterministic, fully synthetic** multi-chromosome cohort + a tiny bundled
FAVOR database, so the `GLOWanalyses` templates run **out-of-the-box** on a fresh
clone with no external data. The outputs are *mechanically* valid but
**statistically meaningless** (real-science validation is a separate effort) — the
example validates the template **mechanics** (every stage emits a schema-valid,
non-empty artifact).

This cohort is **shared-storage-independent**: every example config leaves
`symlinked_shared_root` at its base default (`NULL` = local), so a clone without a
`/project` mount just works.

## What ships (tracked) vs what is generated

**Tracked** (committed; ≈ 0.4 MB total):

| Path | What |
|---|---|
| `make-example-cohort.R` | the deterministic generator (raw inputs) |
| `build-intermediates.sh` | runs the early templates → the pre-built intermediates |
| `check-example.R` | the schema self-test (presence + schema + non-emptiness) |
| `cohort/example_plink/chr{21,22}.{bed,bim,fam}` | raw PLINK (the `00/plink-to-gds` input) |
| `cohort/favor-db/` | the tiny FAVOR DB (`FAVORdatabase_chrsplit.csv` + `chr{21,22}_1.csv`) |
| `cohort/pheno.csv` | `sample_id, outcome(0/1), sex, age, subpop` |
| `cohort/known-snps.csv` | ~20 known SNPs (B training + the PI case subset) |
| `cohort/variant_map.csv` | provenance: every variant's gene/category/MAF |
| `cohort/example_gds_favor/flexible/gds/chr{21,22}_favor.gds` | the pre-built **aGDS** |
| `cohort/b_model.rds`, `cohort/pi_models/` | the pre-built **B model** + **PI ensemble** |

**Generated** (git-ignored; rebuilt by the two scripts — see `.gitignore`):
`cohort/example_gds/`, `cohort/example_pcs/`, `cohort/example_pheno/`,
`cohort/example_ld/`, `cohort/example_gds_favor/flexible/csv/`, and every run's
`runs/example/*/outputs/`.

## Cohort shape

- **2 chromosomes** (21, 22) → exercises per-chr parallelism + per-chr aggregation.
- **300 samples**, **2 sub-populations** (A/B) with slightly shifted allele
  frequencies → non-degenerate PCA + null-model structure; ≈ 50/50 case/control.
- **~370 variants** (~186 chr21, ~188 chr22), **sited inside 8 real chr21/22 genes**
  (`SIK1B, CBSL, U2AF1L5, CRYAA2` on chr21; `POTEH, IL17RA, HDHD5, SLC25A18` on
  chr22, from the bundled `GLOWr` `genes_info`) plus intergenic variants — so
  **gene / window / coding** region definition all resolve with the **default**
  machinery (`GLOWr::define_regions_gene(chr)`). Gene variants carry coding
  categories (plof / missense / disruptive_missense / synonymous / …); the
  intergenic tier supplies the **common-variant** MAF tier (≥ 100 variants at
  MAF ≥ 0.01) that PCA and the single-variant LDSC calibration need.

## Regenerate (one command each)

```bash
conda activate r_env            # run from the GLOWanalyses directory (paths are root-relative)
Rscript data-example/make-example-cohort.R      # raw inputs
bash    data-example/build-intermediates.sh     # pre-built assets
```

`make-example-cohort.R` is deterministic (`set.seed`), so it regenerates the
identical cohort. `build-intermediates.sh` then runs `00/plink-to-gds`,
`00/annotate-favor` (the aGDS + its CSV), `00/compute-pcs`,
`00/assemble-pheno-covar`, and `01-training/{b,pi}`, and copies the B model + PI
ensemble up to `cohort/`.

> **Note on the aGDS.** `00/annotate-favor` (`GLOWr::annotate_favor`) writes the
> `annotation/info/FunctionalAnnotation` node as a **folder of per-feature
> sub-nodes** — the STAARpipeline `gds2agds` format that **both** the per-feature
> scan reads (`seqGetData(..., ".../FunctionalAnnotation/<feature>")`) **and** the
> PI control loader read. The example data-prep config sets `favor_features` to the
> 16 PI features + the 3 string coding nodes
> (`genecode_comprehensive_category`, `genecode_comprehensive_exonic_category`,
> `metasvm_pred`), so the produced aGDS feeds PI training/prediction **and** the
> coding masks directly — no separate aGDS builder is needed.

## Run the demo + self-test

```bash
bash run-example.sh                  # (GLOWanalyses dir) 02 → 03 → 04 over the pre-built assets, then check
REBUILD=1 bash run-example.sh        # also regenerate the cohort + 00/01 intermediates
```

`check-example.R` (also invoked by `run-example.sh`) asserts presence + schema +
non-emptiness for every stage and exits non-zero on any failure. It deliberately
does **not** assert p-value / AUC magnitudes (meaningless on synthetic data).

