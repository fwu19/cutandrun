/*
* Call peaks with both MACS2 and SEACR
*/

include { SEACR_CALLPEAK as SEACR_CALLPEAK_IGG                          } from "../../modules/nf-core/seacr/callpeak/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_NOIGG                        } from "../../modules/nf-core/seacr/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_IGG_NARROW                   } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NOIGG_NARROW                 } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_IGG_BROAD                    } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NOIGG_BROAD                  } from "../../modules/nf-core/macs2/callpeak/main"
include { BEDTOOLS_INTERSECT as SEACR_PEAKS_BEDTOOLS_INTERSECT          } from "../../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT   } from "../../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT    } from "../../modules/nf-core/bedtools/intersect/main"

workflow CALL_PEAKS {
    take:
    bedgraph_markdup
    bedgraph_dedup
    bam_markdup

    main:
    /*
    * MODULE: Call peaks using SEACR with and without IgG control
    */

    /*
     * CHANNEL: Separate bedgraphs into target/control for SEACR
     * for SEACR, use markdup for target and dedup for control
     * CHANNEL: Create target/control pairings
     */
    bedgraph_markdup
        .filter { it -> it[0].is_control == false }
        .map    { it -> [ it[0].control_group, it ] }
        .set { ch_bedgraph_target }
    // ch_bedgraph_target.view()
    // [ control_group, meta, bedgraph ]

    bedgraph_dedup
        .filter { it -> it[0].is_control == true }
        .map { it -> [it[0].control_group, it] }
        .set { ch_bedgraph_control }
    // ch_bedgraph_control.view
    // [ control_group, meta, bedgraph ]

    ch_bedgraph_target
        .filter { it -> it[0] != "" }
        .cross( ch_bedgraph_control )
        .map { it -> [ it[0][1][0], it[0][1][1], it[1][1][1] ] }
        .set { ch_bedgraph_paired }
        // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]

    /*
    * MODULE: Call peaks using SEACR with IgG control
    */
    seacr_peaks                  = Channel.empty()
    seacr_peaks_filtered         = Channel.empty()
    seacr_peaks_igg              = Channel.empty()
    seacr_peaks_noigg            = Channel.empty()

    SEACR_CALLPEAK_IGG (
        ch_bedgraph_paired,
        params.seacr_peak_threshold
    )
    ch_seacr_peaks_igg    = SEACR_CALLPEAK_IGG.out.bed
    ch_versions = SEACR_CALLPEAK_IGG.out.versions
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //SEACR_CALLPEAK_IGG.out.bed | view


    /*
    * CHANNEL: Add fake control channel
    */
    ch_bedgraph_target
        .map{ it -> [ it[1][0], it[1][1], [] ] }
        .set { ch_bedgraph_target_fctrl }
    // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
    // ch_bedgraph_target_fctrl | view

    SEACR_CALLPEAK_NOIGG (
        ch_bedgraph_target_fctrl,
        params.seacr_peak_threshold
    )
    SEACR_CALLPEAK_NOIGG.out.bed
        .branch {
            with_control: it[0].control_group != ""
            no_control: it[0].control_group == ""
        }
        .set { ch_seacr_peaks_noigg}
    ch_versions = ch_versions.mix(SEACR_CALLPEAK_NOIGG.out.versions)
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //SEACR_NO_IGG.out.bed | view


    /*
    * CHANNEL: mix igg and noigg SEACR peaks for cases with control
    */

    ch_seacr_peaks_igg
        .join ( ch_seacr_peaks_noigg.with_control )
        .map { it -> [it[0], it[1], it[2]]}
        .set { ch_seacr_peaks_intersect }

    SEACR_PEAKS_BEDTOOLS_INTERSECT(
        ch_seacr_peaks_intersect,
        [[:],[]]
    )
    ch_versions = ch_versions.mix(SEACR_PEAKS_BEDTOOLS_INTERSECT.out.versions)
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //SEACR_PEAKS_BEDTOOLS_INTERSECT.out.intersect | view


    /*
    * Call MACS2 peaks with and without IgG control
    */

    /*
    * CHANNEL: Separate bams into target/control
    * for MACS2 use markdup for both target and control
    */
    bam_markdup
        .map { it -> [ it[0].control_group, it ] }
        .branch {
            target: it[1][0].is_control == false
            control: it[1][0].is_control == true
        }
        .set { ch_bam }

    //ch_bam.target | view
    // [ control_group, meta, bam ]
    //ch_bam.control | view
    // [ control_group, meta, bam ]


    /*
    * CHANNEL: Create target/control pairings
    */
    ch_bam.target
        .filter { it -> it[0] != "" }
        .cross ( ch_bam.control )
        .map {
            row -> [ row[0][1][0], row[0][1][1], row[1][1][1] ]
        }
        .set { ch_bam_paired }
    //ch_bam_paired | view
    // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BAM, CONTROL_BAM]


    ch_macs2_peaks_narrow_filtered  = Channel.empty()
    ch_macs2_peaks_narrow           = Channel.empty()
    ch_macs2_peaks_igg_narrow       = Channel.empty()
    ch_macs2_peaks_noigg_narrow     = Channel.empty()

    ch_macs2_peaks_broad            = Channel.empty()
    ch_macs2_peaks_broad_filtered   = Channel.empty()
    ch_macs2_peaks_igg_broad        = Channel.empty()
    ch_macs2_peaks_noigg_broad      = Channel.empty()

    MACS2_CALLPEAK_IGG_NARROW (
        ch_bam_paired,
        params.macs_gsize
    )
    ch_macs2_peaks_igg_narrow       = MACS2_CALLPEAK_IGG_NARROW.out.peak
    ch_peaks_summits_igg_narrow     = MACS2_CALLPEAK_IGG_NARROW.out.bed
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_IGG_NARROW.out.versions)

    MACS2_CALLPEAK_IGG_BROAD (
        ch_bam_paired,
        params.macs_gsize
    )
    ch_macs2_peaks_igg_broad       = MACS2_CALLPEAK_IGG_BROAD.out.peak
    ch_peaks_summits_igg_broad     = MACS2_CALLPEAK_IGG_BROAD.out.bed
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_IGG_BROAD.out.versions)

    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //MACS2_CALLPEAK_IGG.out.peak | view



    /*
    * CHANNEL: Add fake control channel
    */
    ch_bam.target.map{ it -> [ it[1][0], it[1][1], [] ] }
    .set { ch_bam_target_fctrl }
    //ch_bam_target_fctrl | view
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]

    MACS2_CALLPEAK_NOIGG_NARROW (
        ch_bam_target_fctrl,
        params.macs_gsize
    )
    MACS2_CALLPEAK_NOIGG_NARROW.out.peak
        .branch {
            with_control: it[0].control_group != ""
            no_control: it[0].control_group == ""
        }
        .set { ch_macs2_peaks_noigg_narrow }
    //ch_peaks_summits_noigg_narrow     = MACS2_CALLPEAK_NOIGG_NARROW.out.bed
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_NOIGG_NARROW.out.versions)

    MACS2_CALLPEAK_NOIGG_BROAD (
        ch_bam_target_fctrl,
        params.macs_gsize
    )
    // MACS2_CALLPEAK_NOIGG.out.peak | view
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    MACS2_CALLPEAK_NOIGG_BROAD.out.peak
        .branch {
            with_control: it[0].control_group != ""
            no_control: it[0].control_group == ""
        }
        .set { ch_macs2_peaks_noigg_broad }
    //ch_peaks_summits_noigg_broad     = MACS2_CALLPEAK_NOIGG_BROAD.out.bed
    ch_versions = ch_versions.mix(MACS2_CALLPEAK_NOIGG_NARROW.out.versions)



    /*
    * CHANNEL: pair igg and noigg MACS2 narrow peaks
    */

    ch_macs2_peaks_igg_narrow
        .join ( ch_macs2_peaks_noigg_narrow.with_control )
        .map { row -> [row[0], row[1], row[2]] }
        .set { ch_macs2_peaks_narrow_intersect }

    MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT(
        ch_macs2_peaks_narrow_intersect,
        [[:],[]]
    )
    ch_versions = ch_versions.mix(MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.versions)
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.intersect | view

    /*
    * CHANNEL: mix igg and noigg MACS2 broad peaks
    */

    ch_macs2_peaks_igg_broad
        .join (ch_macs2_peaks_noigg_broad.with_control )
        .map {row -> [row[0], row[1], row[2]]}
        .set { ch_macs2_peaks_broad_intersect }

    MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT(
        ch_macs2_peaks_broad_intersect,
        [[:],[]]
    )
    ch_versions = ch_versions.mix(MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.versions)
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.intersect | view


    MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.intersect
        .concat(
                ch_macs2_peaks_igg_narrow,
                MACS2_CALLPEAK_NOIGG_NARROW.out.peak,
                MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.intersect,
                ch_macs2_peaks_igg_broad,
                MACS2_CALLPEAK_NOIGG_BROAD.out.peak,
                SEACR_PEAKS_BEDTOOLS_INTERSECT.out.intersect,
                ch_seacr_peaks_igg,
                SEACR_CALLPEAK_NOIGG.out.bed
        )
        .groupTuple(by: 0)
        .set { ch_peaks_all}

    MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.intersect
        .mix    ( ch_macs2_peaks_noigg_narrow.no_control )
        .concat (
            MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.intersect
                .mix( ch_macs2_peaks_noigg_broad.no_control ),
            SEACR_PEAKS_BEDTOOLS_INTERSECT.out.intersect
                .mix( ch_seacr_peaks_noigg.no_control )
        )
        .groupTuple(by: 0)
        .set { ch_peaks_final}

    emit:
    versions = ch_versions
    peaks_all = ch_peaks_all
    peaks_final = ch_peaks_final

}
