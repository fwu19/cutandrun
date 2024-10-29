process ORIGINAL_PEAK_METRICS {
    label "process_single"

    tag "Collect QC of original peaks"

    module = ['fhR/4.1.2-foss-2021b']

    input:
    path( read_metrics )
    path( peak_list, stageAs: "peaks/*" )
    path( rip_list, stageAs: "rip/*" )

    output:
    path( "*.{rds,csv}" ), emit: metrics

    script:
    """
    original_peak_metrics.r

    """

    // input files: read_metrics.csv fragment_length.rds peaks/*.narrowPeak peaks/*.broadPeak peak/*.bed rip/*.reads_in_peak.csv

}
