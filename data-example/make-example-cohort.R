#!/usr/bin/env Rscript
# ============================================================================
# data-example/make-example-cohort.R - deterministic synthetic example cohort
#
# Emits the TRACKED raw inputs of the GLOWanalyses synthetic demo cohort into
# data-example/cohort/ so a fresh clone can run all five stages out-of-the-box.
# Self-contained (reuses the simulation PATTERNS of
# the GLOWr test helper helper-create-test-agds.R - genotype draw, the
# FAVOR annotation-node set, sample annotations - but does NOT source the test
# helper). Deterministic: a fixed set.seed makes the whole cohort regenerable.
#
# OUTPUTS (all under data-example/cohort/, all tracked):
#   example_plink/chr{21,22}.{bed,bim,fam}   per-chr PLINK (the 00/plink-to-gds input)
#   favor-db/FAVORdatabase_chrsplit.csv      FAVOR chunk map (custom, covers chr21/22)
#   favor-db/chr{21,22}_1.csv                FAVOR chunks (variant_vcf-keyed; 16 PI feats
#                                            + coding nodes), the 00/annotate-favor input
#   pheno.csv                                sample_id, outcome (0/1), sex, age, subpop
#   known-snps.csv                           ~20 known SNPs for B + a PI case subset
#   variant_map.csv                          (provenance) every variant's gene/category/MAF
#
# WHY THESE SHAPES:
#   - 2 chromosomes (21, 22) -> per-chr parallelism + per-chr aggregation are
#     genuinely exercised.
#   - ~300 samples, 2 sub-populations w/ slightly different allele freqs -> PCA +
#     the null model have non-degenerate structure.
#   - ~180 variants/chr: ~80% rare (MAF 0.001-0.01) for the SNV-set tests +
#     ~35 common (MAF 0.05-0.30, a few in mild LD) so compute_pcs_gds' prune+PCA
#     is non-degenerate.
#   - variants SITED INSIDE real chr21/22 genes from the bundled GLOWr genes_info
#     so gene / window / coding region definition all resolve with the DEFAULT
#     machinery (GLOWr::define_regions_gene(chr)); a few intergenic for windows.
#
# Usage (from project root):
#   conda activate r_env
#   Rscript data-example/make-example-cohort.R

suppressMessages({
  library(data.table)   # the only external dep; PLINK is hand-written below
})

set.seed(20260608L)  # deterministic: regenerates the identical cohort

# ---- Paths (repo-relative; run from project root) ----
root     <- "data-example/cohort"
plink_dir <- file.path(root, "example_plink")
favor_dir <- file.path(root, "favor-db")
for (d in c(root, plink_dir, favor_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 1. Sample design: ~300 samples, 2 sub-populations, ~50/50 case/control
# ===========================================================================
n_samples <- 300L
sample_ids <- sprintf("ID%04d", seq_len(n_samples))
# Two sub-populations (A, B) with slightly different allele freqs (PCA structure).
subpop <- rep(c("A", "B"), length.out = n_samples)
# Covariates: sex (binary), age (continuous). Outcome ~50/50, mild subpop effect
# so the null model + scan are signal-shaped (not degenerate) but meaningless.
sex <- rbinom(n_samples, 1L, 0.5)
age <- round(rnorm(n_samples, mean = 55, sd = 10), 1)
# Outcome: logit with a small subpop + sex effect, then balanced-ish Bernoulli.
lin <- -0.2 + 0.4 * (subpop == "B") + 0.2 * sex + scale(age)[, 1] * 0.1
outcome <- rbinom(n_samples, 1L, plogis(lin))

pheno <- data.frame(
  sample_id = sample_ids,
  outcome   = as.integer(outcome),
  sex       = as.integer(sex),
  age       = age,
  subpop    = subpop,
  stringsAsFactors = FALSE)

# ===========================================================================
# 2. Gene scaffold: real chr21/22 genes from the bundled GLOWr genes_info
#    (so define_regions_gene(chr) resolves the synthetic variants by default).
# ===========================================================================
# Picked small protein-coding genes (3-40 kb) from GLOWr::genes_info. We hardcode
# the coordinates here so the generator does NOT depend on GLOWr being loadable;
# they are the literal bundled-table values (verified at build time on chr21/22).
gene_scaffold <- rbind(
  data.frame(gene = "SIK1B",    chr = 21L, start = 6111134L, end = 6123739L),
  data.frame(gene = "CBSL",     chr = 21L, start = 6444869L, end = 6468040L),
  data.frame(gene = "U2AF1L5",  chr = 21L, start = 6484623L, end = 6499261L),
  data.frame(gene = "CRYAA2",   chr = 21L, start = 6560714L, end = 6564489L),
  data.frame(gene = "POTEH",    chr = 22L, start = 15690026L, end = 15721631L),
  data.frame(gene = "IL17RA",   chr = 22L, start = 17084954L, end = 17115694L),
  data.frame(gene = "HDHD5",    chr = 22L, start = 17137511L, end = 17165287L),
  data.frame(gene = "SLC25A18", chr = 22L, start = 17563439L, end = 17590994L),
  stringsAsFactors = FALSE)

# Coding-category palette per gene-variant (so coding masks find >=1 set/category
# on the >=2-variant common support). Drawn so plof/missense/synonymous/... all
# appear. Each entry is c(GENCODE.Category, GENCODE.EXONIC.Category, MetaSVM).
.coding_palette <- list(
  stopgain      = c("exonic", "stopgain", ""),
  stoploss      = c("exonic", "stoploss", ""),
  splicing      = c("splicing", "", ""),
  missense_D    = c("exonic", "nonsynonymous SNV", "D"),  # disruptive missense
  missense_T    = c("exonic", "nonsynonymous SNV", "T"),  # tolerated missense
  synonymous    = c("exonic", "synonymous SNV", ""),
  intronic      = c("intronic", "", ""),
  upstream      = c("upstream", "", ""))

# ===========================================================================
# 3. Per-chromosome variant simulation
# ===========================================================================
# The 16 PI features the FAVOR DB must carry (GLOWr:::.default_PI_features()).
pi_features <- c(
  "apc_conservation_v2", "apc_epigenetics", "apc_epigenetics_active",
  "apc_epigenetics_repressed", "apc_epigenetics_transcription",
  "apc_protein_function_v3", "apc_local_nucleotide_diversity_v3",
  "apc_mutation_density", "apc_transcription_factor", "apc_mappability",
  "apc_proximity_to_tsstes", "apc_proximity_to_coding_v2", "apc_micro_rna",
  "cadd_phred", "linsight", "fathmm_xf")

# Build one chromosome's variants: site most inside its genes (coding palette),
# the rest intergenic (for windows). Returns a list of per-variant records.
simulate_chr <- function(chr) {
  genes_chr <- gene_scaffold[gene_scaffold$chr == chr, , drop = FALSE]
  recs <- list()
  vidx <- 0L

  # --- Gene-resident variants: ~28 per gene, coding categories cycled ---
  # cats: 2 plof-ish (stopgain/stoploss/splicing) + several missense + several
  # synonymous + intronic/upstream, so each category has >=2 variants per gene.
  cat_plan <- c(
    rep("stopgain", 2L), rep("stoploss", 1L), rep("splicing", 2L),
    rep("missense_D", 5L), rep("missense_T", 4L),
    rep("synonymous", 6L), rep("intronic", 5L), rep("upstream", 3L))  # 28
  for (gi in seq_len(nrow(genes_chr))) {
    g <- genes_chr[gi, ]
    n_g <- length(cat_plan)
    # Positions spread across the gene body (interior, avoids exact boundaries).
    pos_g <- as.integer(round(seq(g$start + 50L, g$end - 50L, length.out = n_g)))
    for (k in seq_len(n_g)) {
      vidx <- vidx + 1L
      cat_k <- cat_plan[k]
      pal <- .coding_palette[[cat_k]]
      # MAF: gene-resident variants are mostly rare; a few synonymous made common
      # to seed common-tier LD without polluting plof/missense rarity.
      is_common <- (cat_k == "synonymous" && k %% 3L == 0L)
      maf_k <- if (is_common) runif(1, 0.05, 0.20) else runif(1, 0.001, 0.01)
      recs[[length(recs) + 1L]] <- list(
        chr = chr, pos = pos_g[k], gene = g$gene, category = cat_k,
        gencode = pal[1], exonic = pal[2], metasvm = pal[3],
        maf = maf_k, common = is_common)
    }
  }

  # --- Intergenic variants (between/after genes): for the window scan + the
  # common-variant tier the PCA and the single-variant LDSC calibration need.
  # ~75 spread across the chromosome span outside genes, in three MAF tiers:
  #   ~40% COMMON (0.05-0.30), ~25% LOW-FREQ (0.02-0.05), ~35% RARE (0.001-0.01).
  # The common + low-freq tiers give >=100 variants at MAF >= 0.01 across both
  # chromosomes, which is the LDSC regression's minimum (estimate_inflation_factor
  # keeps MAF >= 0.01; ldsc_regression needs >= 100 such pairs).
  span_lo <- min(genes_chr$start) - 200000L
  span_hi <- max(genes_chr$end)   + 200000L
  n_inter <- 80L
  inter_pos <- as.integer(round(seq(span_lo, span_hi, length.out = n_inter + 2L)))
  inter_pos <- inter_pos[-c(1L, length(inter_pos))]  # drop the exact endpoints
  # Drop any that fall inside a gene (keep them strictly intergenic).
  in_gene <- vapply(inter_pos, function(p)
    any(p >= genes_chr$start & p <= genes_chr$end), logical(1))
  inter_pos <- inter_pos[!in_gene]
  for (k in seq_along(inter_pos)) {
    vidx <- vidx + 1L
    tier <- k %% 5L
    if (tier %in% c(0L, 1L)) {           # ~40% common
      maf_k <- runif(1, 0.05, 0.30); is_common <- TRUE
    } else if (tier == 2L) {             # ~20% low-frequency (>= 0.02)
      maf_k <- runif(1, 0.02, 0.05); is_common <- TRUE
    } else {                              # ~40% rare
      maf_k <- runif(1, 0.001, 0.01); is_common <- FALSE
    }
    recs[[length(recs) + 1L]] <- list(
      chr = chr, pos = inter_pos[k], gene = NA_character_, category = "intergenic",
      gencode = "intergenic", exonic = "", metasvm = "",
      maf = maf_k, common = is_common)
  }

  recs
}

records <- c(simulate_chr(21L), simulate_chr(22L))
# De-duplicate any coincident positions within a chromosome (keep first).
key <- vapply(records, function(r) paste(r$chr, r$pos, sep = "-"), character(1))
records <- records[!duplicated(key)]

# ===========================================================================
# 4. Genotypes: binomial(2, MAF) per variant, with mild LD for some commons
#    and a slight per-subpopulation frequency shift (PCA + null-model structure).
# ===========================================================================
n_variants <- length(records)
geno <- matrix(0L, nrow = n_variants, ncol = n_samples)  # variants x samples
is_B <- (subpop == "B")
# Track a "previous common variant per chr" to induce mild LD (copy + noise).
prev_common_geno <- list("21" = NULL, "22" = NULL)
for (v in seq_len(n_variants)) {
  r <- records[[v]]
  # Subpop-specific MAF: shift B by a small delta (bounded to (0, 0.5)).
  delta <- if (r$common) 0.04 else 0.004
  maf_A <- min(max(r$maf - delta / 2, 1e-4), 0.5)
  maf_B <- min(max(r$maf + delta / 2, 1e-4), 0.5)
  g <- integer(n_samples)
  g[!is_B] <- rbinom(sum(!is_B), 2L, maf_A)
  g[is_B]  <- rbinom(sum(is_B),  2L, maf_B)
  # Mild LD: for a common variant, mix in ~40% of the previous common variant's
  # genotype on this chr (a few correlated pairs -> non-degenerate LD prune).
  ck <- as.character(r$chr)
  if (r$common) {
    pg <- prev_common_geno[[ck]]
    if (!is.null(pg)) {
      swap <- which(rbinom(n_samples, 1L, 0.4) == 1L)
      g[swap] <- pg[swap]
    }
    prev_common_geno[[ck]] <- g
  }
  geno[v, ] <- g
}

# Guard: every variant must be polymorphic in the cohort (MAC >= 1) so no stage
# trips an emptiness/monomorphic guard. Re-draw any monomorphic variant with a
# couple of forced minor alleles.
for (v in seq_len(n_variants)) {
  if (sum(geno[v, ]) == 0L) {
    j <- sample.int(n_samples, 2L)
    geno[v, j] <- 1L
  }
}

# ===========================================================================
# 5. Variant metadata frame + FAVOR annotation values
# ===========================================================================
chr_v <- vapply(records, `[[`, integer(1), "chr")
pos_v <- vapply(records, `[[`, integer(1), "pos")
gene_v <- vapply(records, function(r) as.character(r$gene), character(1))
cat_v  <- vapply(records, `[[`, character(1), "category")
gencode_v <- vapply(records, `[[`, character(1), "gencode")
exonic_v  <- vapply(records, `[[`, character(1), "exonic")
metasvm_v <- vapply(records, `[[`, character(1), "metasvm")
maf_target <- vapply(records, `[[`, numeric(1), "maf")
# REF/ALT: simple biallelic SNVs (A/G odd index, C/T even) - single-base, SNV.
ref_v <- ifelse(seq_len(n_variants) %% 2L == 1L, "A", "C")
alt_v <- ifelse(seq_len(n_variants) %% 2L == 1L, "G", "T")
# Realized cohort MAF (for provenance + known-SNP table).
maf_realized <- rowMeans(geno) / 2
maf_realized <- pmin(maf_realized, 1 - maf_realized)

variant_map <- data.frame(
  chr = chr_v, pos = pos_v, ref = ref_v, alt = alt_v,
  gene = gene_v, category = cat_v,
  gencode = gencode_v, exonic = exonic_v, metasvm = metasvm_v,
  maf_target = round(maf_target, 5), maf_realized = round(maf_realized, 5),
  stringsAsFactors = FALSE)
fwrite(variant_map, file.path(root, "variant_map.csv"))

# --- FAVOR feature values: deterministic random draws on plausible scales ---
# (Statistically meaningless; only the SCHEMA + non-NA-ness are validated.)
draw_feat <- function(name) {
  switch(name,
    cadd_phred = runif(n_variants, 0, 40),
    linsight   = runif(n_variants, 0, 1),
    fathmm_xf  = runif(n_variants, 0, 1),
    # All apc_* on a 0-50-ish APC scale.
    runif(n_variants, 0, 50))
}
feat_values <- lapply(pi_features, draw_feat)
names(feat_values) <- pi_features

# ===========================================================================
# 6. Write PLINK per chromosome (hand-rolled PLINK1 binary, SNP-major)
# ===========================================================================
# We write the .bed/.bim/.fam directly rather than via SNPRelate::snpgdsGDS2BED,
# because SeqArray::seqBED2GDS (which GLOWr::plink_to_gds calls) fails on
# snpgdsGDS2BED-produced .bed files in this SeqArray build ("'dimidx' should have 3
# element(s)"). A directly-written canonical PLINK1 .bed reads cleanly through
# seqBED2GDS, so the 00/plink-to-gds template runs unmodified. See
#
# .bed convention: declare A1 (bim col 5) = ALT, A2 (bim col 6) = REF, so PLINK's
# "count of A1" equals the ALT dosage we simulated. PLINK 2-bit codes (per sample):
#   ALT dosage 2 -> 00 (hom A1A1), 1 -> 10 (het), 0 -> 11 (hom A2A2); 01 = missing.
# (SeqArray then reports $ref = REF, $alt = ALT, and $dosage = REF count = 2 - ALT;
# MAF/MAC are direction-invariant, so the genotypes are faithfully preserved.)
.alt_to_plink_code <- function(dose_alt)
  c("2" = 0L, "1" = 2L, "0" = 3L)[as.character(dose_alt)]  # 00, 10, 11

write_plink_chr <- function(chr) {
  sel <- which(chr_v == chr)
  ord <- sel[order(pos_v[sel])]            # PLINK wants position-sorted SNPs
  g_chr <- geno[ord, , drop = FALSE]       # variants x samples (ALT dosage)
  n_snp <- length(ord)
  out_prefix <- file.path(plink_dir, sprintf("chr%d", chr))

  # --- .bed (magic 6c 1b 01 = PLINK1 SNP-major) ---
  bps <- ceiling(n_samples / 4L)           # bytes per SNP row
  bed <- file(paste0(out_prefix, ".bed"), "wb")
  on.exit(close(bed), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), bed)
  for (s in seq_len(n_snp)) {
    codes <- .alt_to_plink_code(g_chr[s, ])
    packed <- integer(bps)
    for (i in seq_len(n_samples)) {
      bi <- (i - 1L) %/% 4L + 1L
      shift <- ((i - 1L) %% 4L) * 2L
      packed[bi] <- packed[bi] + bitwShiftL(codes[i], shift)
    }
    writeBin(as.raw(packed), bed)
  }
  close(bed); on.exit()  # close now; clear the handler

  # --- .bim: chr  rsID  0(cM)  pos  A1(=ALT)  A2(=REF) ---
  rs_ids <- sprintf("rs%d%06d", chr, seq_len(n_snp))
  bim <- data.frame(
    chr = chr, rsid = rs_ids, cm = 0L, pos = pos_v[ord],
    A1 = alt_v[ord], A2 = ref_v[ord], stringsAsFactors = FALSE)
  data.table::fwrite(bim, paste0(out_prefix, ".bim"), sep = "\t", col.names = FALSE)

  # --- .fam: FID IID PID MID SEX PHENO (sample-aligned) ---
  # The pipeline reads outcome/sex from pheno.csv (not the .fam); these encode a
  # well-formed PLINK file for completeness.
  fam <- data.frame(
    FID = "0", IID = sample_ids, PID = 0L, MID = 0L,
    SEX = ifelse(sex == 1L, 1L, 2L),         # 1 = male, 2 = female
    PHENO = ifelse(outcome == 1L, 2L, 1L),   # 1 = control, 2 = case
    stringsAsFactors = FALSE)
  data.table::fwrite(fam, paste0(out_prefix, ".fam"), sep = " ", col.names = FALSE)
  n_snp
}
n_chr21 <- write_plink_chr(21L)
n_chr22 <- write_plink_chr(22L)

# ===========================================================================
# 7. Write the tiny FAVOR database (chrsplit map + per-chr chunk CSVs)
# ===========================================================================
# One chunk per chromosome covering [min_pos, max_pos]; keyed by variant_vcf.
# Columns match GLOWr::annotate_favor's expectation: variant_vcf, chromosome,
# position, ref_vcf, alt_vcf, <16 PI features>, + the coding nodes.
write_favor <- function() {
  split_rows <- list()
  for (chr in c(21L, 22L)) {
    sel <- which(chr_v == chr)
    p   <- pos_v[sel]
    chunk <- data.frame(
      variant_vcf = paste(chr_v[sel], pos_v[sel], ref_v[sel], alt_v[sel], sep = "-"),
      chromosome  = chr_v[sel],
      position    = pos_v[sel],
      ref_vcf     = ref_v[sel],
      alt_vcf     = alt_v[sel],
      stringsAsFactors = FALSE)
    for (f in pi_features) chunk[[f]] <- feat_values[[f]][sel]
    # Coding annotation columns (so a downstream consumer could also build coding
    # masks from the FAVOR CSV if desired; the aGDS carries them as sub-nodes).
    chunk$genecode_comprehensive_category        <- gencode_v[sel]
    chunk$genecode_comprehensive_exonic_category <- exonic_v[sel]
    chunk$metasvm_pred                           <- metasvm_v[sel]
    chunk$rsid <- sprintf("rs%d%06d", chr, seq_along(sel))
    chunk <- chunk[order(chunk$position), , drop = FALSE]
    fwrite(chunk, file.path(favor_dir, sprintf("chr%d_1.csv", chr)))
    split_rows[[length(split_rows) + 1L]] <- data.frame(
      Chr = chr, File_No = 1L, Start_Pos = min(p), End_Pos = max(p),
      stringsAsFactors = FALSE)
  }
  split <- do.call(rbind, split_rows)
  fwrite(split, file.path(favor_dir, "FAVORdatabase_chrsplit.csv"))
}
write_favor()

# ===========================================================================
# 8. Write pheno.csv + known-snps.csv
# ===========================================================================
fwrite(pheno, file.path(root, "pheno.csv"))

# Known-SNP table for B training (+ a PI case subset): ~20 rows drawn from the
# simulated variants (a mix of rare + common, both chromosomes). Columns are the
# standard names prepare_B_training_data / prepare_PI_case_data auto-detect.
known_idx <- sort(c(
  head(which(chr_v == 21L & maf_realized > 0.02), 6L),
  head(which(chr_v == 21L & maf_realized <= 0.02), 4L),
  head(which(chr_v == 22L & maf_realized > 0.02), 6L),
  head(which(chr_v == 22L & maf_realized <= 0.02), 4L)))
known_idx <- known_idx[!is.na(known_idx)]
n_known <- length(known_idx)
known <- data.frame(
  rsID       = sprintf("rsKNOWN%04d", seq_len(n_known)),
  CHR        = chr_v[known_idx],
  POS        = pos_v[known_idx],
  REF        = ref_v[known_idx],
  ALT        = alt_v[known_idx],
  MAF        = round(maf_realized[known_idx], 4),
  # Plausible GWAS-style summary stats (meaningless; schema only).
  P          = signif(runif(n_known, 1e-12, 1e-4), 3),
  N          = rep(3000L, n_known),
  BETA       = round(rnorm(n_known, 0, 0.3), 4),
  TRAIT_TYPE = rep("binary", n_known),
  FIRST_AUTHOR = rep(c("Demo A", "Demo B"), length.out = n_known),
  stringsAsFactors = FALSE)
fwrite(known, file.path(root, "known-snps.csv"))

# ===========================================================================
# 9. Report
# ===========================================================================
cat("=== Synthetic example cohort written ===\n")
cat(sprintf("  samples        : %d (subpop A/B = %d/%d; cases = %d)\n",
            n_samples, sum(subpop == "A"), sum(subpop == "B"), sum(outcome)))
cat(sprintf("  variants       : %d total (chr21 = %d, chr22 = %d)\n",
            n_variants, n_chr21, n_chr22))
cat(sprintf("  common (MAF>=.05): %d ; rare: %d\n",
            sum(maf_realized >= 0.05), sum(maf_realized < 0.05)))
cat(sprintf("  genes (scaffold): %d (%s)\n", nrow(gene_scaffold),
            paste(gene_scaffold$gene, collapse = ", ")))
cat(sprintf("  known SNPs      : %d (B training + PI case subset)\n", n_known))
cat(sprintf("  PLINK           : %s/chr{21,22}.{bed,bim,fam}\n", plink_dir))
cat(sprintf("  FAVOR DB        : %s (chrsplit + chr{21,22}_1.csv; %d PI feats + coding nodes)\n",
            favor_dir, length(pi_features)))
cat(sprintf("  pheno / known   : %s/{pheno.csv, known-snps.csv}\n", root))
cat("Next: bash data-example/build-intermediates.sh\n")
