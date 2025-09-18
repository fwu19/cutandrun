process ORIGINAL_PEAK_WIDTHS {
    label "process_single"

    tag "Collect peak widths from ${meta.id}"

    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    input:
    tuple val(meta), path( "peaks/*" )

    output:
    tuple val(meta), path( "*.csv" ), emit: csv

    script:
    """
    original_peak_widths.r ${meta.id}

    """

    // input files: read_metrics.csv fragment_length.rds peaks/*.narrowPeak peaks/*.broadPeak peak/*.bed rip/*.reads_in_peak.csv
    // output files: [id].original_peaks.csv, [id].original_peak_width.rds
}
