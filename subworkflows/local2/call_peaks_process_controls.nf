/*
* Call peaks with both MACS2 and SEACR
*/

include { SEACR_CALLPEAK as SEACR_CALLPEAK_IGG                          } from "../../modules/nf-core/seacr/callpeak/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_NOIGG                        } from "../../modules/nf-core/seacr/callpeak/main"
include { MACS2_CALLPEAK as MACS2_NARROW_IGG                   } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_NARROW_NOIGG                 } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_BROAD_IGG                    } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_BROAD_NOIGG                  } from "../../modules/nf-core/macs2/callpeak/main"
include { BEDTOOLS_INTERSECT as SEACR_BEDTOOLS_INTERSECT          } from "../../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as NARROW_BEDTOOLS_INTERSECT   } from "../../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as BROAD_BEDTOOLS_INTERSECT    } from "../../modules/nf-core/bedtools/intersect/main"
include { FILTER_PEAKS as NARROW_FILTER                                        } from "../../modules/local2/filter_peaks"
include { FILTER_PEAKS as BROAD_FILTER                                        } from "../../modules/local2/filter_peaks"


workflow CALL_PEAKS_PROCESS_CONTROLS {
    take:
    bedgraph_markdup
    bedgraph_dedup
    bam_markdup
    igg_dir
    use_igg

    main:
    /*
    * MODULE: Call peaks using SEACR with and without IgG control
    */

    ch_bedgraph_paired              = Channel.empty()
    ch_seacr_peaks                  = Channel.empty()
    ch_seacr_peaks_filtered         = Channel.empty()
    ch_seacr_peaks_igg              = Channel.empty()
    ch_seacr_peaks_noigg            = Channel.empty()

    /*
     * CHANNEL: Separate bedgraphs into target/control for SEACR
     * for SEACR, use markdup for target and dedup for control
     * CHANNEL: Create target/control pairings
     */
    bedgraph_markdup
        .filter { it -> it[0].is_control == false }
        .set { ch_bedgraph_target }
    // ch_bedgraph_target.view()
    // [ meta, bedgraph ]


    /*
    * CHANNEL: Add fake control channel and call SEACR peaks with no igg
    */
    ch_bedgraph_target
        .map{ it -> [ it[0], it[1], [] ] }
        .set { ch_bedgraph_target_fctrl }
    // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
    // ch_bedgraph_target_fctrl | view

    SEACR_CALLPEAK_NOIGG (
        ch_bedgraph_target_fctrl,
        params.seacr_peak_threshold
    )
    ch_seacr_peaks_noigg = SEACR_CALLPEAK_NOIGG.out.bed
    ch_versions = SEACR_CALLPEAK_NOIGG.out.versions
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //SEACR_NO_IGG.out.bed | view

    /*
    * CHANNEL: Call SEACR peaks with igg
    */
    if (use_igg){
        ch_bedgraph_target
            .combine(
                Channel.fromPath( "${igg_dir}/IgG.dedup.unnormalized.bedGraph", checkIfExists: true )
            )
            .set { ch_bedgraph_paired }
        // ch_bedgraph_paired.view()
        // EXAMPLE CHANNEL STRUCT: [[TARGET META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]

        SEACR_CALLPEAK_IGG (
            ch_bedgraph_paired,
            params.seacr_peak_threshold
        )
        ch_seacr_peaks_igg    = SEACR_CALLPEAK_IGG.out.bed
        ch_versions = ch_versions.mix(SEACR_CALLPEAK_IGG.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //SEACR_CALLPEAK_IGG.out.bed | view

        /*
        * CHANNEL: mix igg and noigg SEACR peaks for cases with control
        */
        SEACR_BEDTOOLS_INTERSECT(
            ch_seacr_peaks_igg
            .join ( ch_seacr_peaks_noigg.ifEmpty([[:],[]]) )
            .map { it -> [it[0], it[1], it[2]]},
            [[:],[]]
        )
        ch_seacr_peaks_filtered = SEACR_BEDTOOLS_INTERSECT.out.intersect
        ch_versions = ch_versions.mix(SEACR_BEDTOOLS_INTERSECT.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //SEACR_PEAKS_BEDTOOLS_INTERSECT.out.intersect | view

    }




    /*
    * Call MACS2 peaks with and without IgG control
    */

    /*
    * CHANNEL: Separate bams into target/control
    * for MACS2 use markdup for both target and control
    */
    bam_markdup
        .filter { it -> it[0].is_control == false }
        .set { ch_bam_target }

    //ch_bam_target | view
    // [ meta, bam ]

    ch_bam_paired                   = Channel.empty()

    ch_macs2_peaks_narrow_filtered  = Channel.empty()
    ch_macs2_peaks_narrow           = Channel.empty()
    ch_macs2_peaks_igg_narrow       = Channel.empty()
    ch_macs2_peaks_noigg_narrow     = Channel.empty()

    ch_macs2_peaks_broad            = Channel.empty()
    ch_macs2_peaks_broad_filtered   = Channel.empty()
    ch_macs2_peaks_igg_broad        = Channel.empty()
    ch_macs2_peaks_noigg_broad      = Channel.empty()

    /*
    * CHANNEL: Add fake control channel and call MACS2 peaks with no igg
    */
    ch_bam_target.map{ it -> [ it[0], it[1], [] ] }
    .set { ch_bam_target_fctrl }
    //ch_bam_target_fctrl | view
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]

    MACS2_NARROW_NOIGG (
        ch_bam_target_fctrl,
        params.macs_gsize
    )
    ch_macs2_peaks_noigg_narrow = MACS2_NARROW_NOIGG.out.peak
    //ch_peaks_summits_noigg_narrow     = MACS2_NARROW_NOIGG.out.bed
    ch_versions = ch_versions.mix(MACS2_NARROW_NOIGG.out.versions)

    MACS2_BROAD_NOIGG (
        ch_bam_target_fctrl,
        params.macs_gsize
    )
    // MACS2_CALLPEAK_NOIGG.out.peak | view
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    ch_macs2_peaks_noigg_broad = MACS2_BROAD_NOIGG.out.peak
    //ch_peaks_summits_noigg_broad     = MACS2_BROAD_NOIGG.out.bed
    ch_versions = ch_versions.mix(MACS2_BROAD_NOIGG.out.versions)

    /*
    * Call MACS2 peaks with igg
    */
    if (use_igg){
        ch_bam_target
            .combine(
                Channel.fromPath( "${igg_dir}/IgG.target.markdup.sorted.bam" , checkIfExists: true )
            )
            .set { ch_bam_paired }
        //ch_bam_paired | view
        // EXAMPLE CHANNEL STRUCT: [[TARGET META], TARGET_BAM, CONTROL_BAM]

        MACS2_NARROW_IGG (
            ch_bam_paired,
            params.macs_gsize
        )
        ch_macs2_peaks_igg_narrow       = MACS2_NARROW_IGG.out.peak
        ch_peaks_summits_igg_narrow     = MACS2_NARROW_IGG.out.bed
        ch_versions = ch_versions.mix(MACS2_NARROW_IGG.out.versions)

        MACS2_BROAD_IGG (
            ch_bam_paired,
            params.macs_gsize
        )
        ch_macs2_peaks_igg_broad       = MACS2_BROAD_IGG.out.peak
        ch_peaks_summits_igg_broad     = MACS2_BROAD_IGG.out.bed
        ch_versions = ch_versions.mix(MACS2_BROAD_IGG.out.versions)

        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //MACS2_CALLPEAK_IGG.out.peak | view

        /*
        * CHANNEL: pair igg and noigg MACS2 narrow peaks when available
        */
        NARROW_BEDTOOLS_INTERSECT(
            ch_macs2_peaks_igg_narrow
            .join ( ch_macs2_peaks_noigg_narrow.ifEmpty([[:],[]]) )
            .map { row -> [row[0], row[1], row[2]] },
            [[:],[]]
        )
        ch_macs2_peaks_narrow_filtered =  NARROW_BEDTOOLS_INTERSECT.out.intersect
        ch_versions = ch_versions.mix(NARROW_BEDTOOLS_INTERSECT.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.intersect | view

        NARROW_FILTER(
            ch_macs2_peaks_narrow_filtered.ifEmpty([[:],[]])
        )
        ch_macs2_peaks_narrow_filtered =  NARROW_FILTER.out.file
        ch_versions = ch_versions.mix(NARROW_FILTER.out.versions)

        /*
        * CHANNEL: pair igg and noigg MACS2 broad peaks
        */
        BROAD_BEDTOOLS_INTERSECT(
            ch_macs2_peaks_noigg_broad.ifEmpty([[:],[]])
            .join ( ch_macs2_peaks_igg_broad )
            .map {row -> [row[0], row[1], row[2]]},
            [[:],[]]
        )
        ch_macs2_peaks_broad_filtered =  BROAD_BEDTOOLS_INTERSECT.out.intersect
        ch_versions = ch_versions.mix(BROAD_BEDTOOLS_INTERSECT.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.intersect | view

        BROAD_FILTER(
            ch_macs2_peaks_broad_filtered.ifEmpty([[:],[]])
        )
        ch_macs2_peaks_broad_filtered =  BROAD_FILTER.out.file
        ch_versions = ch_versions.mix(BROAD_FILTER.out.versions)

    }

    /*
    * prepare output channel
    */
    if (use_igg){
        ch_macs2_peaks_narrow_filtered
            .concat(
                ch_macs2_peaks_igg_narrow,
                ch_macs2_peaks_noigg_narrow,
                ch_macs2_peaks_broad_filtered,
                ch_macs2_peaks_igg_broad,
                ch_macs2_peaks_noigg_broad,
                ch_seacr_peaks_filtered,
                ch_seacr_peaks_igg,
                ch_seacr_peaks_noigg
            )
            .groupTuple(by: 0)
            .set { ch_peaks_all}
        // ch_peaks_all.view()
        // [ [meta], [path(peak1), path(peak2), ...] ]


        ch_macs2_peaks_narrow_filtered
            .concat (
            ch_macs2_peaks_broad_filtered,
            ch_seacr_peaks_filtered
            )
            .groupTuple(by: 0)
            .set { ch_peaks_final}
        // ch_peaks_final.view()
        // [ [meta], path(macs2_narrow_peak), path(macs2_broad_peak), path(seacr_peak) ]

    }else{
        ch_macs2_peaks_noigg_narrow
            .concat(
                ch_macs2_peaks_noigg_broad,
                ch_seacr_peaks_noigg
            )
            .groupTuple(by: 0)
            .set { ch_peaks_all}
        // ch_peaks_all.view()
        // [ [meta], [path(peak1), path(peak2), ...] ]
        ch_peaks_final = ch_peaks_all

    }

    emit:
    versions = ch_versions
    peaks_all = ch_peaks_all
    peaks_final = ch_peaks_final

}
