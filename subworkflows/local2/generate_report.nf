/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { PLOT_METRICS              } from '../../modules/local2/plot_metrics'

workflow GENERATE_REPORT {
    take:
    read_metrics
    original_peaks
    replicated_peaks
    consensus_peaks

    main:

    /* Make plots for report */
    PLOT_METRICS(
        read_metrics,
        original_peaks,
        replicated_peaks,
        consensus_peaks
    )


    emit:

}
