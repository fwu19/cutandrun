/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { READS_IN_PEAK             } from '../../modules/local2/reads_in_peak'
include { ORIGINAL_PEAKS                } from '../../modules/local2/original_peaks'
include { ORIGINAL_PEAK_WIDTHS          } from '../../modules/local2/original_peak_widths'
include { REPLICATED_PEAKS              } from '../../modules/local2/replicated_peaks'
include { CONSENSUS_PEAKS               } from '../../modules/local2/consensus_peaks'
include { ANNOTATE_CONSENSUS_PEAKS      } from '../../modules/local2/annotate_consensus_peaks'
include { READS_IN_PEAK as READS_IN_PEAK_COMB_IGG             } from '../../modules/local2/reads_in_peak'
include { ORIGINAL_PEAKS as ORIGINAL_PEAKS_COMB_IGG                } from '../../modules/local2/original_peaks'
include { ORIGINAL_PEAK_WIDTHS as ORIGINAL_PEAK_WIDTHS_COMB_IGG          } from '../../modules/local2/original_peak_widths'
include { REPLICATED_PEAKS as REPLICATED_PEAKS_COMB_IGG              } from '../../modules/local2/replicated_peaks'
include { CONSENSUS_PEAKS as CONSENSUS_PEAKS_COMB_IGG               } from '../../modules/local2/consensus_peaks'
include { ANNOTATE_CONSENSUS_PEAKS as ANNOTATE_CONSENSUS_PEAKS_COMB_IGG      } from '../../modules/local2/annotate_consensus_peaks'

workflow QC_PEAKS {
    take:
    peaks_seacr_filtered
    peaks_seacr_igg
    peaks_seacr_noigg
    peaks_seacr_comb_igg_filtered
    peaks_seacr_comb_igg

    peaks_narrow_filtered
    peaks_narrow_igg
    peaks_narrow_noigg
    peaks_narrow_comb_igg_filtered
    peaks_narrow_comb_igg

    peaks_broad_filtered
    peaks_broad_igg
    peaks_broad_noigg
    peaks_broad_comb_igg_filtered
    peaks_broad_comb_igg

    samplesheet
    samplesheet_comb_igg
    ch_samtools_bam
    use_igg
    min_replicates
    fasta
    gtf
    srcdir


    main:
    ch_versions = Channel.empty()
        ch_orig_csv = Channel.empty()
        ch_orig_widths = Channel.empty()
        ch_rip = Channel.empty()
        ch_rep_bed = Channel.empty()
        ch_rep_csv = Channel.empty()
        ch_conp_bed = Channel.empty()
        ch_conp_csv = Channel.empty()
        ch_conp_ann = Channel.empty()
        ch_orig_csv_comb_igg = Channel.empty()
        ch_orig_widths_comb_igg = Channel.empty()
        ch_rip_comb_igg = Channel.empty()
        ch_rep_bed_comb_igg = Channel.empty()
        ch_rep_csv_comb_igg = Channel.empty()
        ch_conp_bed_comb_igg = Channel.empty()
        ch_conp_csv_comb_igg = Channel.empty()
        ch_conp_ann_comb_igg = Channel.empty()

    if ('individual' in use_igg){
        // get pre-generated peak files
        if ( !params.run_peak_calling ) {
            peaks_seacr_igg = Channel.fromPath("${srcdir}/csv/peaks.individual_igg.seacr_igg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.seacr_igg}", checkIfExists: true) ]}
            peaks_seacr_noigg = Channel.fromPath("${srcdir}/csv/peaks.seacr_noigg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.seacr_noigg}", checkIfExists: true) ]}
            peaks_seacr_filtered = Channel.fromPath("${srcdir}/csv/peaks.individual_igg.seacr_filtered.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.seacr_filtered}", checkIfExists: true) ]}
            peaks_narrow_igg = Channel.fromPath("${srcdir}/csv/peaks.individual_igg.macs2_narrow_igg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_narrow_igg}", checkIfExists: true) ]}
            peaks_narrow_noigg = Channel.fromPath("${srcdir}/csv/peaks.macs2_narrow_noigg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_narrow_noigg}", checkIfExists: true) ]}
            peaks_narrow_filtered = Channel.fromPath("${srcdir}/csv/peaks.individual_igg.macs2_narrow_filtered.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_narrow_filtered}", checkIfExists: true) ]}
            peaks_broad_igg = Channel.fromPath("${srcdir}/csv/peaks.individual_igg.macs2_broad_igg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_broad_igg}", checkIfExists: true) ]}
            peaks_broad_noigg = Channel.fromPath("${srcdir}/csv/peaks.macs2_broad_noigg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_broad_noigg}", checkIfExists: true) ]}
            peaks_broad_filtered = Channel.fromPath("${srcdir}/csv/peaks.individual_igg.macs2_broad_filtered.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_broad_filtered}", checkIfExists: true) ]}

        }

        peaks_narrow_filtered
            .concat(
                peaks_narrow_igg,
                peaks_narrow_noigg,
                peaks_broad_filtered,
                peaks_broad_igg,
                peaks_broad_noigg,
                peaks_seacr_filtered,
                peaks_seacr_igg,
                peaks_seacr_noigg
        )
        .groupTuple(by: 0)
        .set { ch_peaks_all}
        // [ [meta], [path(peak1), path(peak2), ...] ]


        peaks_narrow_filtered
            .mix( peaks_narrow_noigg.filter(it -> it[0].control_group == "") )
            .concat (
                peaks_broad_filtered
                    .mix( peaks_broad_noigg.filter(it -> it[0].control_group == "") ),
                peaks_seacr_filtered
                    .mix( peaks_seacr_noigg.filter(it -> it[0].control_group == "") )
            )
            .groupTuple(by: 0)
            .set { ch_peaks_final}
        // [ [meta], path(macs2_narrow_peak), path(macs2_broad_peak), path(seacr_peak) ]

        /*
        * Compute reads in peak
        */
        ch_peaks_final
            .map { it -> [ it[0].id, it[0], it[1] ]}
            .join(
                ch_samtools_bam.map{ it -> [ it[0].id, it[1] ]}
            )
            .map { it -> [ it[1], it[2], it[3] ]}
            .set{ch_peaks_bam}

        READS_IN_PEAK(
            ch_peaks_bam
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
        // path(peaks_metrics)

        /*
        * Collect peak widths for final peaks
        */

        ORIGINAL_PEAK_WIDTHS(
            ch_peaks_final
        )
        ch_orig_widths = ORIGINAL_PEAK_WIDTHS.out.csv
        // ch_orig_widths.view()
        // path(peaks_widths)

        /*
        * Generate replicated peaks and collect metrics
        */
        REPLICATED_PEAKS(
                ch_peaks_final
                    .filter { it[0].call_rep_peak?.toString()?.toBoolean() }
                    .map { it -> [ [it[0].group, it[0].target], it[1] ]}
                    .groupTuple (by: 0)
                    .map { it -> [ it[0][0], it[0][1], it[1].flatten() ] },
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
                .map { it -> [ it[0], it[1].flatten() ] }
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
                fasta,
                gtf,
                ch_conp_bed.map{it[1]}.flatten()
        )
        ch_conp_ann = ANNOTATE_CONSENSUS_PEAKS.out.txt
        // ch_conp_ann.view()
    }

    if ([ 'group', 'all', 'custom' ].any { it in use_igg }){
        if ( !params.run_peak_calling ) {
            peaks_seacr_comb_igg = Channel.fromPath("${srcdir}/csv/peaks.combine_igg.seacr_igg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.seacr_igg}", checkIfExists: true) ]}
            peaks_seacr_noigg = Channel.fromPath("${srcdir}/csv/peaks.seacr_noigg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.seacr_noigg}", checkIfExists: true) ]}
            peaks_seacr_comb_igg_filtered = Channel.fromPath("${srcdir}/csv/peaks.combine_igg.seacr_filtered.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.seacr_filtered}", checkIfExists: true) ]}
            peaks_narrow_comb_igg = Channel.fromPath("${srcdir}/csv/peaks.combine_igg.macs2_narrow_igg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_narrow_igg}", checkIfExists: true) ]}
            peaks_narrow_noigg = Channel.fromPath("${srcdir}/csv/peaks.macs2_narrow_noigg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_narrow_noigg}", checkIfExists: true) ]}
            peaks_narrow_comb_igg_filtered = Channel.fromPath("${srcdir}/csv/peaks.combine_igg.macs2_narrow_filtered.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_narrow_filtered}", checkIfExists: true) ]}
            peaks_broad_comb_igg = Channel.fromPath("${srcdir}/csv/peaks.combine_igg.macs2_broad_igg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_broad_igg}", checkIfExists: true) ]}
            peaks_broad_noigg = Channel.fromPath("${srcdir}/csv/peaks.macs2_broad_noigg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_broad_noigg}", checkIfExists: true) ]}
            peaks_broad_comb_igg_filtered = Channel.fromPath("${srcdir}/csv/peaks.combine_igg.macs2_broad_filtered.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.take(it.size() - 1), file("${srcdir}/${it.macs2_broad_filtered}", checkIfExists: true) ]}

        }

        peaks_narrow_comb_igg_filtered
            .concat(
                peaks_narrow_comb_igg,
                peaks_narrow_noigg,
                peaks_broad_comb_igg_filtered,
                peaks_broad_comb_igg,
                peaks_broad_noigg,
                peaks_seacr_comb_igg_filtered,
                peaks_seacr_comb_igg,
                peaks_seacr_noigg
        )
        .groupTuple(by: 0)
        .set { ch_peaks_all_comb_igg}
        // [ [meta], [path(peak1), path(peak2), ...] ]


        peaks_narrow_comb_igg_filtered
            .mix( peaks_narrow_noigg.filter(it -> it[0].control_group == "") )
            .concat (
                peaks_broad_comb_igg_filtered
                    .mix( peaks_broad_noigg.filter(it -> it[0].control_group == "") ),
                peaks_seacr_comb_igg_filtered
                    .mix( peaks_seacr_noigg.filter(it -> it[0].control_group == "") )
            )
            .groupTuple(by: 0)
            .set { ch_peaks_final_comb_igg}
        // [ [meta], path(macs2_narrow_peak), path(macs2_broad_peak), path(seacr_peak) ]

        /*
        * Compute reads in peak
        */
        ch_peaks_final_comb_igg
            .map { it -> [ it[0].id, it[0], it[1] ]}
            .join(
                ch_samtools_bam.map{ it -> [ it[0].id, it[1] ]}
            )
            .map { it -> [ it[1], it[2], it[3] ]}
            .set{ch_peaks_bam_comb_igg}

        READS_IN_PEAK_COMB_IGG(
            ch_peaks_bam_comb_igg
        )
        ch_rip_comb_igg = READS_IN_PEAK_COMB_IGG.out.csv
        // ch_rip.view()


        /*
        * Collect metrics for original peaks
        */
        ORIGINAL_PEAKS_COMB_IGG(
            ch_peaks_all_comb_igg
        )
        ch_orig_csv_comb_igg = ORIGINAL_PEAKS_COMB_IGG.out.csv
        // ch_orig_peaks.view()
        // path(peaks_metrics)

        /*
        * Collect peak widths for final peaks
        */

        ORIGINAL_PEAK_WIDTHS_COMB_IGG(
            ch_peaks_final_comb_igg
        )
        ch_orig_widths_comb_igg = ORIGINAL_PEAK_WIDTHS_COMB_IGG.out.csv
        // ch_orig_widths.view()
        // path(peaks_widths)

        /*
        * Generate replicated peaks and collect metrics
        */
        REPLICATED_PEAKS_COMB_IGG(
                ch_peaks_final_comb_igg
                    .filter { it[0].call_rep_peak?.toString()?.toBoolean() }
                    .map { it -> [ [it[0].group, it[0].target], it[1] ]}
                    .groupTuple (by: 0)
                    .map { it -> [ it[0][0], it[0][1], it[1].flatten() ] },
                min_replicates
        )
        ch_rep_bed_comb_igg = REPLICATED_PEAKS_COMB_IGG.out.bed
        ch_rep_csv_comb_igg = REPLICATED_PEAKS_COMB_IGG.out.csv
        //ch_rep_bed.view()
        // [ target, [peaks] ]

        /*
        * Generate consensus peaks and collect metrics
        */
        CONSENSUS_PEAKS_COMB_IGG(
                ch_rep_bed_comb_igg
                .groupTuple ( by: 0 )
                .map { it -> [ it[0], it[1].flatten() ] }
                .combine ( samplesheet_comb_igg )
        )
        ch_conp_bed_comb_igg = CONSENSUS_PEAKS_COMB_IGG.out.bed
        ch_conp_csv_comb_igg = CONSENSUS_PEAKS_COMB_IGG.out.csv
        // ch_conp_bed.view()
        // [ target, path(conp) ]

        /*
        * Annotate consensus peaks
        */
        ANNOTATE_CONSENSUS_PEAKS_COMB_IGG(
                fasta,
                gtf,
                ch_conp_bed_comb_igg.map{it[1]}.flatten()
        )
        ch_conp_ann_comb_igg = ANNOTATE_CONSENSUS_PEAKS_COMB_IGG.out.txt
        // ch_conp_ann.view()
    }

    emit:
    rip = ch_rip
    orig_csv = ch_orig_csv
    orig_widths = ch_orig_widths
    rep_bed = ch_rep_bed
    rep_csv = ch_rep_csv
    conp_bed = ch_conp_bed
    conp_csv = ch_conp_csv
    conp_ann = ch_conp_ann
    rip_comb_igg = ch_rip_comb_igg
    orig_csv_comb_igg = ch_orig_csv_comb_igg
    orig_widths_comb_igg = ch_orig_widths_comb_igg
    rep_bed_comb_igg = ch_rep_bed_comb_igg
    rep_csv_comb_igg = ch_rep_csv_comb_igg
    conp_bed_comb_igg = ch_conp_bed_comb_igg
    conp_csv_comb_igg = ch_conp_csv_comb_igg
    conp_ann_comb_igg = ch_conp_ann_comb_igg
    versions = ch_versions
}
