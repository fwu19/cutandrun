/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { READS_IN_PEAK             } from '../../modules/local2/reads_in_peak'
include { ORIGINAL_PEAK_METRICS     } from '../../modules/local2/original_peak_metrics'
include { REPLICATED_PEAKS          } from '../../modules/local2/replicated_peaks'
include { CONSENSUS_PEAKS           } from '../../modules/local2/consensus_peaks'

workflow PEAK_QC {
    take:
    input
    read_metrics
    bam
    peaks_all
    peaks_final

    main:
    /* Compute reads in peak */
    peaks_all
        .join(bam)
        .set{ch_peak_bam}
    READS_IN_PEAK(
        ch_peak_bam
    )

    /*
    * Collect metrics for original peaks
    */
    ORIGINAL_PEAK_METRICS(
        read_metrics,
        peaks_all.map{ it -> it[1] }.flatten().collect(),
        READS_IN_PEAK.out.rip.collect{it[1]}
    )

    // original_peaks.view()


    /*
    * Generate replicated peaks and collect metrics
    */
    REPLICATED_PEAKS(
        input,
        ORIGINAL_PEAK_METRICS.out.data,
        peaks_final.collect{it[1]},
        params.min_replicates
    )


    /*
    * Generate consensus peaks and collect metrics
    */
    CONSENSUS_PEAKS(
        REPLICATED_PEAKS.out.data
    )

    emit:

    //original_peak_metrics = ORIGINAL_PEAK_METRICS.out.data
    // *.{csv,rds}
    //replicated_peak_metrics = REPLICATED_PEAKS.out.data
    // *.{csv,rds}
    //consensus_peak_metrics = CONSENSUS_PEAKS.out.data
    // *.{csv,rds}
    //consensus_peaks = CONSENSUS_PEAKS.out.peak
    // *.bed

}
