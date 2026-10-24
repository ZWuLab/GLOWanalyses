version 1.0

# ==============================================================================
# File: cloud/terra/glow_scan.wdl
# Purpose: GLOW on Terra, and on any engine that runs WDL in Docker. One workflow, GLOWScan,
#   runs the scan stages of GLOWanalyses as thin tasks over the templates shipped inside the
#   image: with region_type = "single_variant" the per-variant scan with its genomic-control
#   calibration (02-single-variant), and with region_type = gene | window | coding | noncoding
#   the SNV-set scan with its aggregation (03-snv-set, 04-summary). Every task writes a
#   config.R from its inputs, with the field names of the GLOWanalyses configuration reference,
#   and calls the template from /opt/GLOWanalyses, so the analysis logic lives in one place.
#   Inputs: per-chromosome annotated GDS files, the phenotype bundle, the PI ensemble as a
#   tarball, and a B model (or a closed-form B function). The SLURM route and this one write
#   the same tables.
#
#   Layout: every task links its inputs under /tmp/glow (gds/chr<N>.gds, models/, pheno/,
#   calibration/) and links /tmp/glow/run to a directory in its working directory, because the
#   Stage-2 templates take their data paths from the Stage-1 snapshot and every machine must
#   see the same paths; the engine delocalizes what lands under run/outputs/.
#
#   Authors: Zheyang Wu (Worcester Polytechnic Institute), for ZWuLab.
#   Derived from the Terra workflows written by Alberto Brusati (Landers lab, UMass Chan
#   Medical School), 2026: the one-workflow shape with the region_type switch, the input names
#   and the resource formula are his.

workflow GLOWScan {

  meta {
    description: "GLOW scan stages: the single-variant scan with calibration, or the gene, window, coding or non-coding SNV-set scan with aggregation, as thin tasks over the GLOWanalyses templates"
    author: "Zheyang Wu (ZWuLab); derived from the Terra workflows written by Alberto Brusati (Landers lab, UMass Chan Medical School), 2026"
  }

  parameter_meta {
    region_type: "single_variant | gene | window | coding | noncoding. single_variant runs 02-single-variant (scan, LD scores, calibration); the others run 03-snv-set and 04-summary (scan, aggregation, plots)"
    run_label: "A free label; it only enters the returned run_tag (region_type or region_type_label), not the output file names"
    gds_files: "Per-chromosome annotated GDS files (the aGDS of 00-data-prep), parallel to chromosomes; any basenames, the tasks link each one as chr<N>.gds"
    chromosomes: "Chromosome numbers as strings, 1 to 22, distinct, parallel to gds_files (the i-th number names the i-th file)"
    pheno_bundle: "The phenotype bundle .rds: a list with sample_id, Y (the trait), X (the covariate matrix) and trait (binary or continuous), as 00-data-prep/assemble-pheno-covar.R writes it"
    pi_models_tar: "Set-based runs only, required there: the PI ensemble as a .tar.gz of model_*.rds files, at the archive root or inside one directory"
    b_model_rds: "A trained B model .rds (01-training/train-b-model.R); when absent, b_func is used"
    b_func: "A closed-form B function in R (a string), used when b_model_rds is absent"
    rare_maf_cutoff: "Variant filter (filter_spec): keep variants with MAF below this value"
    variant_type: "Variant filter: SNV (the only value coding runs accept), indel or variant"
    min_mac: "Variant filter: minimum minor-allele count"
    min_variants: "Variant filter: minimum number of variants a set needs to be tested"
    categories: "coding runs: the coding categories (plof, plof_ds, missense, disruptive_missense, synonymous, ptv, ptv_ds)"
    noncoding_categories: "noncoding runs: the gene-centric masks (upstream, downstream, UTR, promoter_CAGE, promoter_DHS, enhancer_CAGE, enhancer_DHS, ncRNA); the aGDS must carry the five FAVOR fields the masks read"
    window_size: "window runs: the sliding window length in base pairs"
    step_size: "window runs: the step between windows in base pairs"
    chunk_size_mb: "window runs: the genome is split into chunks of this many megabases; one chunk is one unit of work inside a chromosome's task"
    gene_max_records: "gene runs: a gene with more GDS records in its span than this is tiled into segments (see the GLOWanalyses README, Cores and large genes)"
    gene_segment_bp: "gene runs: the segment length in base pairs for such a gene"
    staar_enabled: "Run the embedded STAAR comparator on GLOW's own variant sets (STAAR_O, ACAT_O and the family omnibi as extra columns; STAAR_O_glowG in coding aggregates). Coding runs with it need a binary trait. The native STAARpipeline stage of the SLURM route is not part of this workflow"
    use_spa: "GLOW's single-variant statistic inside the set-based tests: auto (SPA for a binary trait, standard otherwise), spa, or standard"
    primary_test: "The test the aggregation ranks and thresholds on (GLOW_Omni by default)"
    alpha: "Family-wise level of the aggregation's Bonferroni threshold"
    top_k: "gene runs: rows of glow_top_hits.csv"
    output_format: "csv; fst is accepted for window runs only (the gene, coding and noncoding tables are read as CSV)"
    single_variant_scan: "Calibration of a set-based run: marginal_all.csv of a single_variant run on the same phenotype bundle and covariates. With ld_scores it estimates the inflation factor by calibration_method; without any calibration input the run is uncalibrated (factor 1)"
    ld_scores: "ld_scores.rds of that single_variant run; required with single_variant_scan when calibration_method is ldsc_intercept"
    calibration_method: "ldsc_intercept (the LD score regression intercept; needs ld_scores and at least two chromosomes with enough variants) or lambda_gc (the genomic-control lambda; works with one chromosome). Also the method of a single_variant run's own calibration"
    inflation_factor: "A manual inflation factor for a set-based run instead of the two files above; giving both is an error"
    use_SPA: "single_variant runs: the saddlepoint approximation for a binary trait"
    chunk_size: "single_variant runs: variants per scan chunk"
    mac_cutoff: "single_variant runs: minimum minor-allele count of a tested variant"
    ld_window: "single_variant runs: the LD-score window in base pairs on each side"
    ld_segment: "single_variant runs: core variants per LD-score segment (the BLAS block)"
    ld_sample_n: "single_variant runs: individuals subsampled for the LD scores; 0 = all"
    ld_ref_maf_min: "single_variant runs: only SNVs with MAF at least this value enter the LD scores, as targets and as neighbors; 0 = every SNV. Restricting them changes the estimated intercept, so validate the calibration on the cohort"
    docker_image: "The GLOW image; the default is the one published with this GLOWanalyses version"
    cpu: "Cores of the scan tasks, passed to the templates as n_cores (forked workers per chromosome)"
    mem_gb_prepare: "Memory of the Prepare task (Stage 1 of a set-based run), in GB"
    mem_gb_scan: "Memory of each per-chromosome scan task, in GB; the LD-score step of a single_variant run holds matrices that grow with the sample count and the variants per window"
    mem_gb_aggregate: "Memory of the Calibrate and Aggregate tasks, in GB"
    disk_multiplier: "A task's disk is its localized inputs times this factor, plus disk_extra_gb"
    disk_extra_gb: "Disk margin of every task, in GB"
    preemptible: "Preemptible attempts per task on Terra (0 = none)"
  }

  input {
    String        region_type   = "gene"
    String        run_label     = ""
    Array[File]   gds_files
    Array[String] chromosomes
    File          pheno_bundle
    File?         pi_models_tar
    File?         b_model_rds
    String        b_func        = "function(maf) sqrt(maf * (1 - maf))"

    # The variant filter (filter_spec of base-config/snv_set_base.R)
    Float         rare_maf_cutoff = 0.5
    String        variant_type    = "SNV"
    Int           min_mac         = 1
    Int           min_variants    = 1

    # Set-based settings (the field names of the configuration reference)
    Array[String] categories = ["plof", "plof_ds", "missense", "disruptive_missense", "synonymous"]
    Array[String] noncoding_categories = ["upstream", "downstream", "UTR", "promoter_CAGE", "promoter_DHS", "enhancer_CAGE", "enhancer_DHS", "ncRNA"]
    Int           window_size      = 100000
    Int           step_size        = 10000
    Int           chunk_size_mb    = 5
    Int           gene_max_records = 10000
    Int           gene_segment_bp  = 10000
    Boolean       staar_enabled    = false
    String        use_spa          = "auto"
    String        primary_test     = "GLOW_Omni"
    Float         alpha            = 0.05
    Int           top_k            = 100
    String        output_format    = "csv"

    # Calibration of a set-based run
    File?         single_variant_scan
    File?         ld_scores
    String        calibration_method = "ldsc_intercept"
    Float?        inflation_factor

    # Single-variant settings (base-config/single_variant_base.R)
    Boolean       use_SPA        = true
    Int           chunk_size     = 2000
    Int           mac_cutoff     = 1
    Int           ld_window      = 1000000
    Int           ld_segment     = 2000
    Int           ld_sample_n    = 0
    Float         ld_ref_maf_min = 0.0

    # Resources
    String        docker_image     = "ghcr.io/zwulab/glowanalyses:v0.3.0"
    Int           cpu              = 8
    Int           mem_gb_prepare   = 32
    Int           mem_gb_scan      = 32
    Int           mem_gb_aggregate = 16
    Float         disk_multiplier  = 1.5
    Int           disk_extra_gb    = 20
    Int           preemptible      = 0
  }

  String  tag = if run_label == "" then region_type else region_type + "_" + run_label
  Boolean is_single_variant = region_type == "single_variant"

  if (is_single_variant) {
    scatter (i in range(length(gds_files))) {
      call ScanVariants {
        input:
          gds_file = gds_files[i], chrom = chromosomes[i], pheno_bundle = pheno_bundle,
          use_SPA = use_SPA, chunk_size = chunk_size, mac_cutoff = mac_cutoff,
          ld_window = ld_window, ld_segment = ld_segment, ld_sample_n = ld_sample_n,
          ld_ref_maf_min = ld_ref_maf_min, calibration_method = calibration_method,
          docker_image = docker_image, cpu = cpu, mem_gb = mem_gb_scan,
          disk_multiplier = disk_multiplier, disk_extra_gb = disk_extra_gb, preemptible = preemptible
      }
    }
    call Calibrate {
      input:
        chromosomes = chromosomes, pheno_bundle = pheno_bundle,
        scan_tables = ScanVariants.scan_table, ld_tables = ScanVariants.ld_table,
        ld_window = ld_window, ld_segment = ld_segment, ld_sample_n = ld_sample_n,
        ld_ref_maf_min = ld_ref_maf_min, calibration_method = calibration_method,
        docker_image = docker_image, mem_gb = mem_gb_aggregate,
        disk_multiplier = disk_multiplier, disk_extra_gb = disk_extra_gb, preemptible = preemptible
    }
  }

  if (!is_single_variant) {
    call Prepare {
      input:
        region_type = region_type, gds_files = gds_files, chromosomes = chromosomes,
        pheno_bundle = pheno_bundle, pi_models_tar = pi_models_tar,
        b_model_rds = b_model_rds, b_func = b_func,
        rare_maf_cutoff = rare_maf_cutoff, variant_type = variant_type, min_mac = min_mac,
        min_variants = min_variants, categories = categories,
        noncoding_categories = noncoding_categories, window_size = window_size,
        step_size = step_size, chunk_size_mb = chunk_size_mb, gene_max_records = gene_max_records,
        gene_segment_bp = gene_segment_bp, staar_enabled = staar_enabled, use_spa = use_spa,
        primary_test = primary_test, alpha = alpha, top_k = top_k, output_format = output_format,
        single_variant_scan = single_variant_scan, ld_scores = ld_scores,
        calibration_method = calibration_method, inflation_factor = inflation_factor, cpu = cpu,
        docker_image = docker_image, mem_gb = mem_gb_prepare,
        disk_multiplier = disk_multiplier, disk_extra_gb = disk_extra_gb, preemptible = preemptible
    }
    scatter (i in range(length(gds_files))) {
      call ScanChromosome {
        input:
          region_type = region_type, gds_file = gds_files[i], chrom = chromosomes[i],
          shared = Prepare.shared, pi_models_tar = pi_models_tar,
          b_model_rds = b_model_rds, pheno_bundle = pheno_bundle,
          docker_image = docker_image, cpu = cpu, mem_gb = mem_gb_scan,
          disk_multiplier = disk_multiplier, disk_extra_gb = disk_extra_gb, preemptible = preemptible
      }
    }
    call Aggregate {
      input:
        region_type = region_type, shared = Prepare.shared,
        result_files = flatten(ScanChromosome.results),
        docker_image = docker_image, mem_gb = mem_gb_aggregate,
        disk_multiplier = disk_multiplier, disk_extra_gb = disk_extra_gb, preemptible = preemptible
    }
  }

  output {
    # single_variant
    File?        marginal_all                = Calibrate.marginal_all
    File?        ld_scores_out               = Calibrate.ld_scores_out
    File?        single_variant_calibrated   = Calibrate.calibrated
    File?        single_variant_calibration  = Calibrate.calibration_txt
    File?        single_variant_plots        = Calibrate.plots
    Array[File]? single_variant_logs         = ScanVariants.log
    # set-based
    File?                shared_artifacts    = Prepare.shared
    File?                prepare_log         = Prepare.log
    Array[Array[File]]?  scan_results        = ScanChromosome.results
    Array[Array[File]]?  scan_logs           = ScanChromosome.logs
    Array[File]?         aggregated          = Aggregate.aggregated
    File?                aggregated_tar      = Aggregate.aggregated_tar
    File?                aggregate_log       = Aggregate.log
    String               run_tag             = tag
  }
}

# ------------------------------------------------------------------------------
# Single-variant run: one task per chromosome (marginal-scan.R --chr, compute-ld-scores.R
# --chr), then Calibrate combines the chromosomes and calibrates.
# ------------------------------------------------------------------------------
task ScanVariants {
  input {
    File    gds_file
    String  chrom
    File    pheno_bundle
    Boolean use_SPA
    Int     chunk_size
    Int     mac_cutoff
    Int     ld_window
    Int     ld_segment
    Int     ld_sample_n
    Float   ld_ref_maf_min
    String  calibration_method
    String  docker_image
    Int     cpu
    Int     mem_gb
    Float   disk_multiplier
    Int     disk_extra_gb
    Int     preemptible
  }
  Int disk_gb = ceil((size(gds_file, "GB") + size(pheno_bundle, "GB")) * disk_multiplier) + disk_extra_gb

  command <<<
    set -euo pipefail
    # Input validation (before any work).
    [[ "~{chrom}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || { echo "ERROR: chromosome '~{chrom}' is not one of 1..22" >&2; exit 2; }
    case "~{calibration_method}" in ldsc_intercept|lambda_gc) ;; *) echo "ERROR: calibration_method must be ldsc_intercept or lambda_gc (got '~{calibration_method}')" >&2; exit 2 ;; esac
    export OPENBLAS_NUM_THREADS=1
    TASK=$(pwd)
    mkdir -p "$TASK/run" /tmp/glow/gds /tmp/glow/pheno
    ln -sfn "$TASK/run" /tmp/glow/run
    ln -sfn "~{gds_file}" "/tmp/glow/gds/chr~{chrom}.gds"
    ln -sfn "~{pheno_bundle}" /tmp/glow/pheno/pheno_covar.rds

    cat > /tmp/glow/run/config.R <<'RCONF'
    # config.R written by glow_scan.wdl (ScanVariants); field names as in base-config/single_variant_base.R
    source("base-config/single_variant_base.R", local = TRUE)
    chroms             <- c(~{chrom})
    gds_dir            <- "/tmp/glow/gds"
    gds_pattern        <- "chr{chr}.gds"
    pheno_path         <- "/tmp/glow/pheno/pheno_covar.rds"
    ld_scores_path     <- "/tmp/glow/run/outputs/ld/ld_scores.rds"   # --chr writes ld_chr<N>.rds beside it
    use_SPA            <- ~{if use_SPA then "TRUE" else "FALSE"}
    chunk_size         <- ~{chunk_size}L
    mac_cutoff         <- ~{mac_cutoff}L
    n_cores            <- ~{cpu}L
    ld_window          <- ~{ld_window}L
    ld_segment         <- ~{ld_segment}L
    ld_sample_n        <- ~{ld_sample_n}L
    ld_ref_maf_min     <- ~{ld_ref_maf_min}
    calibration_method <- "~{calibration_method}"
    symlinked_shared_root <- NULL
    RCONF

    cd /opt/GLOWanalyses
    cat /opt/glow-image-manifest.txt
    Rscript 02-single-variant/marginal-scan.R     --config /tmp/glow/run/config.R --chr ~{chrom}
    Rscript 02-single-variant/compute-ld-scores.R --config /tmp/glow/run/config.R --chr ~{chrom}
  >>>

  output {
    File scan_table = "run/outputs/results/marginal_chr~{chrom}.csv"
    File ld_table   = "run/outputs/ld/ld_chr~{chrom}.rds"
    File log        = "run/outputs/logs/scan_log.txt"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
  }
}

task Calibrate {
  input {
    Array[String] chromosomes
    File          pheno_bundle
    Array[File]   scan_tables
    Array[File]   ld_tables
    Int           ld_window
    Int           ld_segment
    Int           ld_sample_n
    Float         ld_ref_maf_min
    String        calibration_method
    String        docker_image
    Int           mem_gb
    Float         disk_multiplier
    Int           disk_extra_gb
    Int           preemptible
  }
  Int disk_gb = ceil((size(scan_tables, "GB") + size(ld_tables, "GB") + size(pheno_bundle, "GB")) * disk_multiplier) + disk_extra_gb

  command <<<
    set -euo pipefail
    # Input validation (before any work).
    CHR=( ~{sep=" " chromosomes} )
    [ "${#CHR[@]}" -gt 0 ] || { echo "ERROR: chromosomes is empty" >&2; exit 2; }
    [ "$(printf '%s\n' "${CHR[@]}" | sort -u | wc -l)" -eq "${#CHR[@]}" ] || { echo "ERROR: chromosomes repeats a value" >&2; exit 2; }
    case "~{calibration_method}" in ldsc_intercept|lambda_gc) ;; *) echo "ERROR: calibration_method must be ldsc_intercept or lambda_gc (got '~{calibration_method}')" >&2; exit 2 ;; esac
    export OPENBLAS_NUM_THREADS=1
    TASK=$(pwd)
    mkdir -p "$TASK/run/outputs/results" "$TASK/run/outputs/ld" /tmp/glow/gds /tmp/glow/pheno
    ln -sfn "$TASK/run" /tmp/glow/run
    ln -sfn "~{pheno_bundle}" /tmp/glow/pheno/pheno_covar.rds
    # The per-chromosome tables are linked, not copied (the templates read them in place).
    for f in ~{sep=" " scan_tables}; do ln -s "$f" "$TASK/run/outputs/results/"; done
    for f in ~{sep=" " ld_tables};   do ln -s "$f" "$TASK/run/outputs/ld/"; done

    cat > /tmp/glow/run/config.R <<'RCONF'
    # config.R written by glow_scan.wdl (Calibrate); field names as in base-config/single_variant_base.R
    source("base-config/single_variant_base.R", local = TRUE)
    chroms             <- c(~{sep=", " chromosomes})
    gds_dir            <- "/tmp/glow/gds"
    gds_pattern        <- "chr{chr}.gds"
    pheno_path         <- "/tmp/glow/pheno/pheno_covar.rds"
    ld_scores_path     <- "/tmp/glow/run/outputs/ld/ld_scores.rds"
    ld_window          <- ~{ld_window}L
    ld_segment         <- ~{ld_segment}L
    ld_sample_n        <- ~{ld_sample_n}L
    ld_ref_maf_min     <- ~{ld_ref_maf_min}
    calibration_method <- "~{calibration_method}"
    symlinked_shared_root <- NULL
    RCONF

    cd /opt/GLOWanalyses
    cat /opt/glow-image-manifest.txt
    Rscript 02-single-variant/marginal-scan.R     --config /tmp/glow/run/config.R --combine
    Rscript 02-single-variant/compute-ld-scores.R --config /tmp/glow/run/config.R --combine
    Rscript 02-single-variant/calibrate.R         --config /tmp/glow/run/config.R --method ~{calibration_method}
  >>>

  output {
    File marginal_all    = "run/outputs/results/marginal_all.csv"
    File ld_scores_out   = "run/outputs/ld/ld_scores.rds"
    File calibrated      = "run/outputs/results/single_variant_all_calibrated.csv"
    File calibration_txt = "run/outputs/results/single_variant_calibration.txt"
    File plots           = "run/outputs/plots/marginal_plots.pdf"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
  }
}

# ------------------------------------------------------------------------------
# Set-based run: Prepare (Stage 1, once), ScanChromosome (Stage 2, one per chromosome),
# Aggregate (Stage 3, once).
# ------------------------------------------------------------------------------
task Prepare {
  input {
    String        region_type
    Array[File]   gds_files
    Array[String] chromosomes
    File          pheno_bundle
    File?         pi_models_tar
    File?         b_model_rds
    String        b_func
    Float         rare_maf_cutoff
    String        variant_type
    Int           min_mac
    Int           min_variants
    Array[String] categories
    Array[String] noncoding_categories
    Int           window_size
    Int           step_size
    Int           chunk_size_mb
    Int           gene_max_records
    Int           gene_segment_bp
    Boolean       staar_enabled
    String        use_spa
    String        primary_test
    Float         alpha
    Int           top_k
    String        output_format
    File?         single_variant_scan
    File?         ld_scores
    String        calibration_method
    Float?        inflation_factor
    Int           cpu
    String        docker_image
    Int           mem_gb
    Float         disk_multiplier
    Int           disk_extra_gb
    Int           preemptible
  }
  Int disk_gb = ceil((size(gds_files, "GB") + size(pheno_bundle, "GB") + size(pi_models_tar, "GB") + size(b_model_rds, "GB")
                      + size(single_variant_scan, "GB") + size(ld_scores, "GB")) * disk_multiplier) + disk_extra_gb

  command <<<
    set -euo pipefail
    # Input validation (before any work): parallel, non-empty, distinct chromosomes 1..22.
    GDS=( ~{sep=" " gds_files} ); CHR=( ~{sep=" " chromosomes} )
    [ "${#CHR[@]}" -gt 0 ] || { echo "ERROR: chromosomes is empty" >&2; exit 2; }
    [ "${#CHR[@]}" -eq "${#GDS[@]}" ] || { echo "ERROR: gds_files (${#GDS[@]}) and chromosomes (${#CHR[@]}) differ in length" >&2; exit 2; }
    for c in "${CHR[@]}"; do [[ "$c" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || { echo "ERROR: chromosome '$c' is not one of 1..22" >&2; exit 2; }; done
    [ "$(printf '%s\n' "${CHR[@]}" | sort -u | wc -l)" -eq "${#CHR[@]}" ] || { echo "ERROR: chromosomes repeats a value" >&2; exit 2; }
    PI="~{default='' pi_models_tar}"
    [ -n "$PI" ] || { echo "ERROR: a set-based run (region_type = ~{region_type}) needs pi_models_tar, the PI ensemble as a .tar.gz of model_*.rds files" >&2; exit 2; }
    case "~{use_spa}" in auto|spa|standard) ;; *) echo "ERROR: use_spa must be auto, spa or standard (got '~{use_spa}')" >&2; exit 2 ;; esac
    case "~{calibration_method}" in ldsc_intercept|lambda_gc) ;; *) echo "ERROR: calibration_method must be ldsc_intercept or lambda_gc (got '~{calibration_method}')" >&2; exit 2 ;; esac
    ~{if defined(single_variant_scan) && defined(inflation_factor) then "echo 'ERROR: give either single_variant_scan (with ld_scores) or inflation_factor, not both' >&2; exit 2" else "true"}
    ~{if defined(single_variant_scan) && !defined(ld_scores) && calibration_method == "ldsc_intercept" then "echo 'ERROR: calibration_method ldsc_intercept needs ld_scores beside single_variant_scan; give ld_scores or choose lambda_gc' >&2; exit 2" else "true"}
    export OPENBLAS_NUM_THREADS=1
    TASK=$(pwd)
    mkdir -p "$TASK/run" /tmp/glow/gds /tmp/glow/models/pi_models /tmp/glow/pheno /tmp/glow/calibration
    ln -sfn "$TASK/run" /tmp/glow/run
    for i in "${!GDS[@]}"; do ln -sfn "${GDS[$i]}" "/tmp/glow/gds/chr${CHR[$i]}.gds"; done
    ln -sfn "~{pheno_bundle}" /tmp/glow/pheno/pheno_covar.rds
    tar -xzf "$PI" -C /tmp/glow/models/pi_models
    find /tmp/glow/models/pi_models -mindepth 2 -name '*.rds' -exec mv -t /tmp/glow/models/pi_models/ {} +
    echo "PI models: $(ls /tmp/glow/models/pi_models/*.rds | wc -l) files"
    ~{if defined(b_model_rds) then "ln -sfn '" + select_first([b_model_rds]) + "' /tmp/glow/models/b_model.rds" else "true"}
    ~{if defined(single_variant_scan) then "ln -sfn '" + select_first([single_variant_scan]) + "' /tmp/glow/calibration/marginal_all.csv" else "true"}
    ~{if defined(ld_scores) then "ln -sfn '" + select_first([ld_scores]) + "' /tmp/glow/calibration/ld_scores.rds" else "true"}

    cat > /tmp/glow/run/config.R <<'RCONF'
    # config.R written by glow_scan.wdl (Prepare); field names as in base-config/snv_set_base.R
    source("base-config/snv_set_base.R", local = TRUE)
    region_type  <- "~{region_type}"
    chroms       <- c(~{sep=", " chromosomes})
    gds_dir      <- "/tmp/glow/gds"
    gds_pattern  <- "chr{chr}.gds"
    pheno_path   <- "/tmp/glow/pheno/pheno_covar.rds"
    pi_model_dir <- "/tmp/glow/models/pi_models"
    ~{if defined(b_model_rds) then 'b_func <- NULL; b_model_path <- "/tmp/glow/models/b_model.rds"' else "b_func <- " + b_func + "; b_model_path <- NULL"}
    filter_spec  <- list(rare_maf_cutoff = ~{rare_maf_cutoff}, variant_type = "~{variant_type}",
                         min_mac = ~{min_mac}L, min_variants = ~{min_variants}L)
    n_cores      <- ~{cpu}L
    use_spa      <- ~{if use_spa == "auto" then "NULL" else if use_spa == "spa" then "TRUE" else "FALSE"}
    staar_enabled <- ~{if staar_enabled then "TRUE" else "FALSE"}
    staar_modes   <- character(0)   # this workflow runs no native STAARpipeline stage; the comparator is staar_enabled
    categories   <- c("~{sep='", "' categories}")
    noncoding_categories <- c("~{sep='", "' noncoding_categories}")
    window_size  <- ~{window_size}L
    step_size    <- ~{step_size}L
    chunk_size_mb <- ~{chunk_size_mb}
    gene_max_records <- ~{gene_max_records}L
    gene_segment_bp  <- ~{gene_segment_bp}L
    output_format     <- "~{output_format}"
    primary_test      <- "~{primary_test}"
    primary_test_glow <- "~{primary_test}"
    alpha        <- ~{alpha}
    top_k        <- ~{top_k}L
    ~{if defined(single_variant_scan) then 'calibration <- list(enabled = TRUE, method = "' + calibration_method + '", single_variant_scan = "/tmp/glow/calibration/marginal_all.csv", ld_scores = ' + (if defined(ld_scores) then '"/tmp/glow/calibration/ld_scores.rds"' else "NULL") + ")" else "# no calibration files given"}
    ~{if defined(inflation_factor) then "inflation_factor <- " + select_first([inflation_factor]) else "# no manual inflation_factor"}
    symlinked_shared_root <- NULL
    RCONF

    cd /opt/GLOWanalyses
    cat /opt/glow-image-manifest.txt
    Rscript 03-snv-set/prepare.R --config /tmp/glow/run/config.R 2>&1 | tee "$TASK/prepare.log"
    tar -czf "$TASK/shared.tar.gz" -C "$TASK" run/config.R run/outputs/shared
    ls -la "$TASK/run/outputs/shared"
  >>>

  output {
    File shared = "shared.tar.gz"
    File log    = "prepare.log"
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
  }
}

task ScanChromosome {
  input {
    String region_type
    File   gds_file
    String chrom
    File   shared
    File?  pi_models_tar
    File?  b_model_rds
    File   pheno_bundle
    String docker_image
    Int    cpu
    Int    mem_gb
    Float  disk_multiplier
    Int    disk_extra_gb
    Int    preemptible
  }
  Int disk_gb = ceil((size(gds_file, "GB") + size(shared, "GB") + size(pi_models_tar, "GB") + size(pheno_bundle, "GB")
                      + size(b_model_rds, "GB")) * disk_multiplier) + disk_extra_gb

  command <<<
    set -euo pipefail
    [[ "~{chrom}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]] || { echo "ERROR: chromosome '~{chrom}' is not one of 1..22" >&2; exit 2; }
    PI="~{default='' pi_models_tar}"
    [ -n "$PI" ] || { echo "ERROR: a set-based run needs pi_models_tar" >&2; exit 2; }
    export OPENBLAS_NUM_THREADS=1
    TASK=$(pwd)
    mkdir -p "$TASK/run" /tmp/glow/gds /tmp/glow/models/pi_models /tmp/glow/pheno
    ln -sfn "$TASK/run" /tmp/glow/run
    tar -xzf "~{shared}" -C "$TASK"
    ln -sfn "~{gds_file}" "/tmp/glow/gds/chr~{chrom}.gds"
    ln -sfn "~{pheno_bundle}" /tmp/glow/pheno/pheno_covar.rds
    tar -xzf "$PI" -C /tmp/glow/models/pi_models
    find /tmp/glow/models/pi_models -mindepth 2 -name '*.rds' -exec mv -t /tmp/glow/models/pi_models/ {} +
    ~{if defined(b_model_rds) then "ln -sfn '" + select_first([b_model_rds]) + "' /tmp/glow/models/b_model.rds" else "true"}

    cd /opt/GLOWanalyses
    cat /opt/glow-image-manifest.txt
    case "~{region_type}" in
      gene)      Rscript 03-snv-set/run-gene.R      --config /tmp/glow/run/config.R --chr ~{chrom} ;;
      coding)    Rscript 03-snv-set/run-coding.R    --config /tmp/glow/run/config.R --chr ~{chrom} ;;
      noncoding) Rscript 03-snv-set/run-noncoding.R --config /tmp/glow/run/config.R --chr ~{chrom} ;;
      window)
        CHUNKS=$(Rscript -e 'gc <- readRDS("/tmp/glow/run/outputs/shared/01-genome_chunks.rds"); cat(gc$chunk_id[as.character(gc$chr) == "~{chrom}"])')
        echo "window chunks of chr~{chrom}: ${CHUNKS}"
        for k in ${CHUNKS}; do
          Rscript 03-snv-set/run-window.R --config /tmp/glow/run/config.R --chunk "$k"
        done ;;
      *) echo "unknown region_type ~{region_type}" >&2; exit 2 ;;
    esac
    ls -la "$TASK"/run/outputs/results*/
  >>>

  output {
    Array[File] results = glob("run/outputs/results*/*")
    Array[File] logs    = glob("run/outputs/logs/log_*")
  }

  runtime {
    docker: docker_image
    cpu: cpu
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
  }
}

task Aggregate {
  input {
    String      region_type
    File        shared
    Array[File] result_files
    String      docker_image
    Int         mem_gb
    Float       disk_multiplier
    Int         disk_extra_gb
    Int         preemptible
  }
  Int disk_gb = ceil((size(shared, "GB") + size(result_files, "GB")) * disk_multiplier) + disk_extra_gb

  command <<<
    set -euo pipefail
    export OPENBLAS_NUM_THREADS=1
    TASK=$(pwd)
    mkdir -p "$TASK/run" /tmp/glow
    ln -sfn "$TASK/run" /tmp/glow/run
    tar -xzf "~{shared}" -C "$TASK"
    case "~{region_type}" in
      gene|window) DEST="$TASK/run/outputs/results" ;;
      *)           DEST="$TASK/run/outputs/results-glow" ;;
    esac
    mkdir -p "$DEST"
    for f in ~{sep=" " result_files}; do ln -s "$f" "$DEST/"; done   # linked, not copied: the aggregator reads them in place
    echo "gathered $(ls "$DEST" | wc -l) result files into $DEST"

    cd /opt/GLOWanalyses
    cat /opt/glow-image-manifest.txt
    Rscript 04-summary/aggregate.R --config /tmp/glow/run/config.R 2>&1 | tee "$TASK/aggregate.log"
    Rscript 04-summary/plots.R     --config /tmp/glow/run/config.R 2>&1 | tee -a "$TASK/aggregate.log"
    tar -czf "$TASK/aggregated.tar.gz" -C "$TASK" run/outputs/aggregated
  >>>

  output {
    Array[File] aggregated     = glob("run/outputs/aggregated/*")
    File        aggregated_tar = "aggregated.tar.gz"
    File        log            = "aggregate.log"
  }

  runtime {
    docker: docker_image
    cpu: 2
    memory: "~{mem_gb} GB"
    disks: "local-disk ~{disk_gb} SSD"
    preemptible: preemptible
  }
}
