/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { ORIGINAL_PEAK_METRICS     } from '../../modules/local2/original_peak_metrics'
include { REPLICATED_PEAKS          } from '../../modules/local2/replicated_peaks'
include { CONSENSUS_PEAKS           } from '../../modules/local2/consensus_peaks'

workflow PEAK_QC {
    take:
    read_metrics
    reads_in_peak
    peaks_all
    peaks_final

    main:

    /* Collect metrics for original peaks */
    ORIGINAL_PEAK_METRICS(
        read_metrics,
        peaks_all.map{ it -> it[1:]  }.collect().flatten(),
        reads_in_peak
    )
    // original_peaks.view()


    /* Generate replicated peaks and collect metrics  */
    REPLICATED_PEAKS(
        params.input,
        peaks_final.map{ it -> it[1] }.collect(),
        params.min_replicates
    )


    /* Generate consensus peaks and collect metrics  */
    CONSENSUS_PEAKS(
        REPLICATED_PEAKS.out.metrics
    )

    emit:
    original_peak_metrics = ORIGINAL_PEAK_METRICS.out.metrics
    // *.{csv,rds}
    replicated_peak_metrics = REPLICATED_PEAKS.out.metrics
    // *.{csv,rds}
    consensus_peak_metrics = CONSENSUS_PEAKS.out.metrics
    // *.{csv,rds}
    consensus_peaks = CONSENSUS_PEAKS.out.peak
    // *.bed

}
