process PLOT_METRICS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Make plots of read and peak metrics "

    input:
    path( read_metrics )
    path( original_peaks)
    path( replicated_peaks )
    path( consensus_peaks )

    output:
    tuple path( "*" )

    script:
    """
    plot_metrics.r
    
    """
}
