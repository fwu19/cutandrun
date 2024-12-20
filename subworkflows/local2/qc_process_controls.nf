/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { READS_IN_PEAK             } from '../../modules/local2/reads_in_peak'
include { ORIGINAL_PEAKS                } from '../../modules/local2/original_peaks'
include { ORIGINAL_PEAK_WIDTHS          } from '../../modules/local2/original_peak_widths'
include { REPLICATED_HICONF_PEAKS              } from '../../modules/local2/replicated_high_confidence_peaks'

workflow QC_PROCESS_CONTROLS {
    take:
    samplesheet
    genome
    ch_samtools_bam
    ch_peaks_all
    ch_peaks_final
    ch_hiconf_peaks

    main:
        ch_orig_csv = Channel.empty()
        ch_orig_widths = Channel.empty()
        ch_rip = Channel.empty()
        ch_rep_bed = Channel.empty()
        ch_rep_csv = Channel.empty()

        /*
        * Compute reads in peak
        */
        ch_peaks_final
            .join(ch_samtools_bam)
            .set{ch_peak_bam}

        READS_IN_PEAK(
            ch_peak_bam
        )
        ch_rip = READS_IN_PEAK.out.csv
        // ch_rip.view()


        /*
        * Collect metrics for original peaks
        */
        ORIGINAL_PEAKS(
            ch_peaks_all
        )
        ch_orig_csv = ORIGINAL_PEAKS.out.csv
        // ch_orig_peaks.view()
        // path(peak_metrics)

        /*
        * Collect peak widths for final peaks
        */

        ORIGINAL_PEAK_WIDTHS(
            ch_peaks_final
        )
        ch_orig_widths = ORIGINAL_PEAK_WIDTHS.out.csv
        // ch_orig_widths.view()
        // path(peak_widths)

        /*
        * Generate replicated peaks and collect metrics
        */
        REPLICATED_HICONF_PEAKS(
                ch_peaks_final
                    .combine(ch_hiconf_peaks)
        )
        ch_rep_csv = REPLICATED_HICONF_PEAKS.out.csv
        //ch_rep_bed.view()
        // [ target, [peaks] ]


    emit:
    rip = ch_rip
    orig_csv = ch_orig_csv
    orig_widths = ch_orig_widths
    rep_csv = ch_rep_csv

}
