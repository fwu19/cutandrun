process ORIGINAL_PEAKS {
    label "process_single"

    tag "Collect metrics of original peaks from ${meta.id}"

    module = ['fhR/4.1.2-foss-2021b']

    input:
    tuple val(meta), path( "peaks/*" )

    output:
    tuple val(meta), path( "*.csv" ), emit: csv

    script:
    """
    original_peaks.r ${meta.id}

    """

    // input files: read_metrics.csv fragment_length.rds peaks/*.narrowPeak peaks/*.broadPeak peak/*.bed rip/*.reads_in_peak.csv
    // output files: [id].original_peaks.csv
}
