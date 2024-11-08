process CONSENSUS_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Generate consensus peaks from $target"

    input:
    path (samplesheet)
    tuple val (target), path ( "peaks/*" )

    output:
    tuple val(target), path ( "*.csv" ), emit: csv
    tuple val(target), path ( "*.rds" ), emit: rds
    tuple val(target), path ( "*.bed" ), emit: bed

    script:
    """
    consensus_peaks.r $samplesheet $target

    """

    // input files: replicated_peaks.rds replicated_peak_metrics.csv
}
