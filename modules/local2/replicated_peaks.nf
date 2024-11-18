process REPLICATED_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Generate replicated peaks from $group"

    input:
    tuple val( group ), val( target ), path ( "peaks/*" )
    val ( min_reps )


    output:
    tuple val(target), path ( "*.csv" ), emit: csv, optional: true
    tuple val(target), path ( "*.rds" ), emit: rds, optional: true
    tuple val(target), path ( "{multiple_replicates,single_replicate}/*.bed" ), emit: bed, optional: true

    script:
    """
    replicated_peaks.r $group $min_reps

    """

    // input files: original_peaks.rds original_peak_metrics.csv
}
