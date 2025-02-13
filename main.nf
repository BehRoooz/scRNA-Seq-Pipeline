#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = "raw_data"
params.output = "results"

process QUALITY_CONTROL {
    input:
    path input_data

    output:
    path "filtered_data"

    script:
    """
    python ${projectDir}/quality_control.py ${input_data} filtered_data
    """
}

process DIMENSIONALITY_REDUCTION {
    input:
    path filtered_data

    output:
    path "reduced_data"

    script:
    """
    python ${projectDir}/dimensionality_reduction.py ${filtered_data} clustered_data
    """
}

process CLUSTERING {
    input:
    path clustered_data

    output:
    path "clustered_data"

    script:
    """
    python ${projectDir}/clustering.py ${clustered_data} clustered_data
    """
}

process CELL_MARKERS {
    input:
    path clustered_data

    output:
    path "cell_markers.csv"

    script:
    """
    python ${projectDir}/cell_markers.py ${clustered_data} cell_markers.csv
    """
}

workflow {
    input_ch = channel.fromPath(params.input)

    qc_out = QUALITY_CONTROL(input_ch)
    dr_out = DIMENSIONALITY_REDUCTION(qc_out)
    clustering_out = CLUSTERING(dr_out)
    CELL_MARKERS(clustering_out)
}

