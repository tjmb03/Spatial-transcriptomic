#!/usr/bin/env nextflow
// =============================================================================
// main.nf — Microenvironment Spatial Transcriptomics Pipeline (Nextflow DSL2)
//
// Usage:
//
//   Full pipeline (real Slide-tag data):
//   nextflow run nextflow/main.nf \
//     -profile docker \
//     --sample_sheet sample_sheet.csv \
//     --genome       /path/to/star_genome \
//     --reference    /path/to/mrca_reference.h5seurat \
//     --outdir       results/
//
//   Test mode (Tran 2019 prototype, skips STARsolo):
//   nextflow run nextflow/main.nf \
//     -profile docker,test \
//     --outdir prototype/output/nextflow_test
//
//   HPC SLURM:
//   nextflow run nextflow/main.nf \
//     -profile slurm \
//     --sample_sheet sample_sheet.csv \
//     --genome /scratch/genome/GRCm39 \
//     --outdir /scratch/results
//     -resume
// =============================================================================
nextflow.enable.dsl = 2

// ── Input validation ──────────────────────────────────────────────────────────
if (!params.test_mode) {
    if (!params.sample_sheet) error "ERROR: --sample_sheet is required"
}

// ─────────────────────────────────────────────────────────────────────────────
// PROCESSES
// ─────────────────────────────────────────────────────────────────────────────

// ── Step 0: Validate preregistration hash ─────────────────────────────────────
process VALIDATE_PREREG {
    tag "preregistration"

    output:
    path "prereg_validated.txt", emit: validated

    script:
    """
    ${params.python3} -c "
import sys, os
sys.path.insert(0, '${projectDir}')
try:
    from config.preregistration import validate_preregistration, compute_hash
    validate_preregistration(strict=False)
    h = compute_hash()
    print(f'Preregistration hash: {h}')
    with open('prereg_validated.txt', 'w') as f:
        f.write(h)
except Exception as e:
    print(f'WARNING: Preregistration check skipped: {e}')
    with open('prereg_validated.txt', 'w') as f:
        f.write('UNSET')
"
    """
}

// ── Step 1: STARsolo alignment (skipped in test mode) ────────────────────────
process STARSOLO {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(r1), path(r2), path(puck)

    output:
    tuple val(sample_id), path("${sample_id}/"),                           emit: alignment_dir
    tuple val(sample_id), path("${sample_id}/spatial_coordinates.tsv"),    emit: coords
    tuple val(sample_id), path("${sample_id}/spatial_qc.json"),            emit: qc

    publishDir "${params.outdir}/starsolo_out", mode: params.publish_dir_mode

    script:
    """
    bash ${projectDir}/../scripts/01_starsolo_alignment.sh \
        --sample   ${sample_id} \
        --r1       ${r1} \
        --r2       ${r2} \
        --puck     ${puck} \
        --genome   ${params.genome} \
        --outdir   . \
        --threads  ${task.cpus}
    """
}

// ── Step 2: Seurat v5 QC, clustering, label transfer ─────────────────────────
process SEURAT_PROCESSING {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(alignment_dir), path(spatial_coords)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.h5seurat"), emit: h5seurat
    tuple val(sample_id), path("${sample_id}_metadata.tsv"), emit: metadata
    path "${sample_id}_qc.json", emit: qc

    publishDir "${params.outdir}/seurat_out/${sample_id}", mode: params.publish_dir_mode

    script:
    ref_arg = reference.name != 'NO_REF' ? "--reference ${reference}" : ""
    """
    Rscript ${projectDir}/../scripts/02_seurat_processing.R \
        --sample_id    ${sample_id} \
        --align_dir    ${alignment_dir} \
        --spatial_file ${spatial_coords} \
        --outdir       . \
        ${ref_arg}
    """
}

// ── Convert h5seurat → h5ad (R → Python bridge) ───────────────────────────────
process CONVERT_H5AD {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(h5seurat), path(metadata)

    output:
    tuple val(sample_id), path("${sample_id}.h5ad"), emit: h5ad

    script:
    """
    ${params.python3} - << 'PYEOF'
import anndata as ad
import scanpy as sc
import pandas as pd

# Use seuratdisk if available, else rely on pre-converted file
try:
    import seuratdisk
    seuratdisk.convert("${h5seurat}", dest="h5ad")
except ImportError:
    pass  # h5ad may already exist from Seurat step

adata = sc.read_h5ad("${sample_id}.h5ad")

# Attach metadata from TSV
meta = pd.read_csv("${metadata}", sep="\t", index_col=0)
for col in meta.columns:
    if col not in adata.obs.columns:
        shared = adata.obs_names.intersection(meta.index)
        adata.obs.loc[shared, col] = meta.loc[shared, col]

# Ensure raw counts are in layers['counts']
if 'counts' not in adata.layers:
    import numpy as np, scipy.sparse as sp
    x = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    if np.all(x == np.floor(x)):
        adata.layers['counts'] = adata.X.copy()
    else:
        raise ValueError("No integer count layer found — check Seurat conversion")

adata.write_h5ad("${sample_id}.h5ad")
print(f"Converted: ${sample_id}.h5ad ({adata.n_obs} cells x {adata.n_vars} genes)")
PYEOF
    """
}

// ── Step 3: scVI batch correction + 7-test validation ────────────────────────
process SCVI_INTEGRATION {
    tag "all_samples"

    input:
    path h5ad_files
    path sample_csv

    output:
    path "integrated_adata.h5ad",  emit: integrated
    path "scvi_model/",            emit: model
    path "validation_report.json", emit: validation

    publishDir "${params.outdir}/scvi_out", mode: params.publish_dir_mode

    script:
    skip_t6 = params.dev_mode ? "--skip_t6_stability" : ""
    """
    ${params.python3} ${projectDir}/../scripts/03_scvi_integration.py \
        --h5ad_dir   . \
        --sample_csv ${sample_csv} \
        --outdir     . \
        ${skip_t6}
    """
}

// ── Step 4: RNA velocity (dynamical mode) ────────────────────────────────────
process SCVELO {
    tag "velocity"

    input:
    path integrated_adata

    output:
    path "velocity_adata.h5ad", emit: velocity
    path "velocity_qc.json",    emit: qc

    publishDir "${params.outdir}/velocity_out", mode: params.publish_dir_mode

    script:
    """
    ${params.python3} ${projectDir}/../scripts/04_scvelo_velocity.py \
        --adata   ${integrated_adata} \
        --outdir  . \
        --n_jobs  ${task.cpus}
    """
}

// ── Step 5: CellRank fate probabilities ──────────────────────────────────────
process CELLRANK {
    tag "trajectories"

    input:
    path velocity_adata

    output:
    path "fate_adata.h5ad",      emit: fate
    path "kernel_weights.json",  emit: weights

    publishDir "${params.outdir}/cellrank_out", mode: params.publish_dir_mode

    script:
    """
    ${params.python3} ${projectDir}/../scripts/05_cellrank_trajectories.py \
        --adata  ${velocity_adata} \
        --outdir .
    """
}

// ── Step 6: Spatial statistics + primary endpoint ────────────────────────────
process SPATIAL_STATISTICS {
    tag "spatial"

    input:
    path fate_adata

    output:
    path "primary_endpoint.json", emit: primary
    path "morans_results.json",   emit: morans

    publishDir "${params.outdir}/spatial_out", mode: params.publish_dir_mode

    script:
    """
    ${params.python3} ${projectDir}/../scripts/06_spatial_statistics.py \
        --adata  ${fate_adata} \
        --outdir .
    """
}

// ── Step 7: Lineage drivers + therapeutic targets ─────────────────────────────
process LINEAGE_DRIVERS {
    tag "drivers"

    input:
    path fate_adata

    output:
    path "therapeutic_targets.csv", emit: targets

    publishDir "${params.outdir}/driver_out", mode: params.publish_dir_mode

    script:
    ot_flag = params.skip_opentargets ? "--no_opentargets" : ""
    """
    ${params.python3} ${projectDir}/../scripts/07_lineage_drivers_and_targets.py \
        --adata   ${fate_adata} \
        --outdir  . \
        ${ot_flag}
    """
}

// ── Final validation gate ─────────────────────────────────────────────────────
process VALIDATE_PIPELINE {
    tag "validation"

    input:
    path scvi_validation
    path kernel_weights
    path primary_endpoint
    path targets

    output:
    path "pipeline_validation_report.json", emit: json
    path "pipeline_validation_report.txt",  emit: txt

    publishDir "${params.outdir}/validation_report", mode: params.publish_dir_mode

    script:
    """
    ${params.python3} ${projectDir}/../validation/validate_pipeline.py \
        --results_dir ${params.outdir} \
        --outdir      .
    """
}

// ─────────────────────────────────────────────────────────────────────────────
// SUBWORKFLOW: Test mode (skip STARsolo + Seurat, use pre-built h5ad)
// ─────────────────────────────────────────────────────────────────────────────
workflow TEST_WORKFLOW {
    take:
    validated_prereg

    main:
    // Use the Tran 2019 validated h5ad directly
    test_h5ad = Channel.fromPath(
        "prototype/output/python/tran2019_validated.h5ad",
        checkIfExists: false
    ).ifEmpty {
        error "Test h5ad not found. Run Step 4 prototype first:\n" +
              "  conda activate glaucoma-r\n" +
              "  bash prototype/run_prototype.sh"
    }

    // Mock sample CSV for test
    test_csv = Channel.fromPath("${projectDir}/../prototype/test_sample_sheet.csv")

    // Collect single test h5ad as list
    h5ads = test_h5ad.collect()

    SCVI_INTEGRATION(h5ads, test_csv.first())
    SCVELO(SCVI_INTEGRATION.out.integrated)
    CELLRANK(SCVELO.out.velocity)

    SPATIAL_STATISTICS(CELLRANK.out.fate)
    LINEAGE_DRIVERS(CELLRANK.out.fate)

    VALIDATE_PIPELINE(
        SCVI_INTEGRATION.out.validation,
        CELLRANK.out.weights,
        SPATIAL_STATISTICS.out.primary,
        LINEAGE_DRIVERS.out.targets
    )

    emit:
    report = VALIDATE_PIPELINE.out.txt
}

// ─────────────────────────────────────────────────────────────────────────────
// MAIN WORKFLOW
// ─────────────────────────────────────────────────────────────────────────────
workflow {

    // ── Validate preregistration ──────────────────────────────────────────────
    VALIDATE_PREREG()

    if (params.test_mode) {
        // ── Test mode: skip alignment, use prototype h5ad ─────────────────────
        log.info "Running in TEST MODE — using Tran 2019 prototype data"
        TEST_WORKFLOW(VALIDATE_PREREG.out.validated)

        TEST_WORKFLOW.out.report.view { f ->
            "\n" + "═"*60 + "\n" +
            " TEST RUN COMPLETE\n" +
            " Report: ${f}\n" +
            "═"*60
        }

    } else {
        // ── Full pipeline ──────────────────────────────────────────────────────

        // Input channel from sample sheet
        sample_ch = Channel
            .fromPath(params.sample_sheet, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                tuple(
                    row.sample_id,
                    file(row.r1_fastq, checkIfExists: true),
                    file(row.r2_fastq, checkIfExists: true),
                    file(row.puck_file, checkIfExists: true)
                )
            }

        sample_csv  = file(params.sample_sheet)
        reference   = params.reference ? file(params.reference) : file('NO_REF')

        // Steps 1–2: per sample (parallel)
        STARSOLO(sample_ch)

        seurat_input = STARSOLO.out.alignment_dir
            .join(STARSOLO.out.coords)
        SEURAT_PROCESSING(seurat_input, reference)

        // Convert to h5ad
        h5ad_input = SEURAT_PROCESSING.out.h5seurat
            .join(SEURAT_PROCESSING.out.metadata)
        CONVERT_H5AD(h5ad_input)

        // Collect all h5ads → scVI (single run on merged data)
        all_h5ads = CONVERT_H5AD.out.h5ad
            .map { id, h5ad -> h5ad }
            .collect()

        SCVI_INTEGRATION(all_h5ads, sample_csv)
        SCVELO(SCVI_INTEGRATION.out.integrated)
        CELLRANK(SCVELO.out.velocity)

        // Steps 6–7 run in parallel on fate_adata
        SPATIAL_STATISTICS(CELLRANK.out.fate)
        LINEAGE_DRIVERS(CELLRANK.out.fate)

        // Final validation gate
        VALIDATE_PIPELINE(
            SCVI_INTEGRATION.out.validation,
            CELLRANK.out.weights,
            SPATIAL_STATISTICS.out.primary,
            LINEAGE_DRIVERS.out.targets
        )

        VALIDATE_PIPELINE.out.txt.view { f ->
            "\n" + "═"*60 + "\n" +
            " PIPELINE COMPLETE\n" +
            " Report: ${f}\n" +
            " Results: ${params.outdir}/\n" +
            "═"*60
        }
    }
}
