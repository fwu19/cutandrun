process CONSENSUS_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Generate consensus peaks and collect metrics "

    input:
    path( replicated_peaks )

    output:
    tuple path( "*" )

    script:
    """
    consensus_peaks.r # replicated_peaks.rds replicated_peak_metrics.csv
    
    """
}
