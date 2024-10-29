process REPLICATED_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Generate replicated peaks and collect metrics"

    input:
    path ( input )
    path ( original_peaks )
    val ( min_reps )


    output:
    path ( "*.{csv,rds}" ), emit: metrics
    path ( "*.bed" ), emit: peak

    script:
    """
    replicated_peaks.r $input $min_reps

    """

    // input files: original_peaks.rds original_peak_metrics.csv
}
