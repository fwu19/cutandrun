/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { READS_IN_PEAK             } from '../../modules/local2/reads_in_peak'
include { ORIGINAL_PEAKS                } from '../../modules/local2/original_peaks'
include { ORIGINAL_PEAK_WIDTHS          } from '../../modules/local2/original_peak_widths'
include { REPLICATED_PEAKS              } from '../../modules/local2/replicated_peaks'
include { CONSENSUS_PEAKS               } from '../../modules/local2/consensus_peaks'
include { ANNOTATE_CONSENSUS_PEAKS      } from '../../modules/local2/annotate_consensus_peaks'

workflow QC_PEAKS {
    take:
    samplesheet
    min_replicates
    genome
    gtf_ann
    ch_samtools_bam
    ch_peaks_all
    ch_peaks_final


    main:
        ch_orig_csv = Channel.empty()
        ch_orig_widths = Channel.empty()
        ch_rip = Channel.empty()
        ch_rep_bed = Channel.empty()
        ch_rep_csv = Channel.empty()
        ch_conp_bed = Channel.empty()
        ch_conp_csv = Channel.empty()
        ch_conp_ann = Channel.empty()
        ch_versions = Channel.empty()

        /*
        * Compute reads in peak
        */
        ch_peaks_final
            .map { it -> [ it[0].id, it[0], it[1] ]}
            .join(
                ch_samtools_bam.map{ it -> [ it[0].id, it[1] ]}
            )
            .map { it -> [ it[1], it[2], it[3] ]}
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
        REPLICATED_PEAKS(
                ch_peaks_final
                    .filter { it[0].call_rep_peak == true }
                    .map { it -> [ [it[0].group, it[0].target], it[1] ]}
                    .groupTuple (by: 0)
                    .map { it -> [ it[0][0], it[0][1], it[1].flatten().collect() ] },
                min_replicates
        )
        ch_rep_bed = REPLICATED_PEAKS.out.bed
        ch_rep_csv = REPLICATED_PEAKS.out.csv
        //ch_rep_bed.view()
        // [ target, [peaks] ]

        /*
        * Generate consensus peaks and collect metrics
        */
        CONSENSUS_PEAKS(
                ch_rep_bed
                .groupTuple ( by: 0 )
                .map { it -> [ it[0], it[1].flatten().collect() ] }
                .combine ( samplesheet )
        )
        ch_conp_bed = CONSENSUS_PEAKS.out.bed
        ch_conp_csv = CONSENSUS_PEAKS.out.csv
        // ch_conp_bed.view()
        // [ target, path(conp) ]

        /*
        * Annotate consensus peaks
        */
        ANNOTATE_CONSENSUS_PEAKS(
                genome,
                gtf_ann,
                ch_conp_bed.collect{it[1]}.flatten()
        )
        ch_conp_ann = ANNOTATE_CONSENSUS_PEAKS.out.txt
        // ch_conp_ann.view()



    emit:
    rip = ch_rip
    orig_csv = ch_orig_csv
    orig_widths = ch_orig_widths
    rep_bed = ch_rep_bed
    rep_csv = ch_rep_csv
    conp_bed = ch_conp_bed
    conp_csv = ch_conp_csv
    conp_ann = ch_conp_ann
    versions = ch_versions
}
