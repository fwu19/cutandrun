/*
* Call peaks with both MACS2 and SEACR
*/

include { MACS2_CALLPEAK as MACS2_NARROW_IGG                   } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_NARROW_NOIGG                 } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_BROAD_IGG                    } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_BROAD_NOIGG                  } from "../../modules/nf-core/macs2/callpeak/main"
include { BEDTOOLS_INTERSECT as NARROW_BEDTOOLS_INTERSECT   } from "../../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as BROAD_BEDTOOLS_INTERSECT    } from "../../modules/nf-core/bedtools/intersect/main"
include { FILTER_PEAKS as NARROW_FILTER                     } from "../../modules/local2/filter_peaks"
include { FILTER_PEAKS as BROAD_FILTER                      } from "../../modules/local2/filter_peaks"

include { MERGE_BAM as MERGE_BAM_TARGET_MARKDUP                      } from "../../modules/local2/merge_bam"
include { MACS2_CALLPEAK as MACS2_NARROW_COMB_IGG                    } from "../../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_BROAD_COMB_IGG                     } from "../../modules/nf-core/macs2/callpeak/main"
include { BEDTOOLS_INTERSECT as NARROW_COMB_IGG_BEDTOOLS_INTERSECT   } from "../../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as BROAD_COMB_IGG_BEDTOOLS_INTERSECT    } from "../../modules/nf-core/bedtools/intersect/main"
include { FILTER_PEAKS as NARROW_COMB_IGG_FILTER                     } from "../../modules/local2/filter_peaks"
include { FILTER_PEAKS as BROAD_COMB_IGG_FILTER                      } from "../../modules/local2/filter_peaks"


workflow CALL_MACS2_PEAKS {
    take:
    bam_markdup
    combine_igg
    meta_combine_igg


    main:

    ch_narrow_igg       = Channel.empty()
    ch_narrow_noigg     = Channel.empty()
    ch_narrow_filtered  = Channel.empty()
    ch_narrow_comb_igg  = Channel.empty()
    ch_narrow_comb_igg_filtered  = Channel.empty()

    ch_broad_igg        = Channel.empty()
    ch_broad_noigg      = Channel.empty()
    ch_broad_filtered   = Channel.empty()
    ch_broad_comb_igg  = Channel.empty()
    ch_broad_comb_igg_filtered  = Channel.empty()


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

    ch_bam.control
        .cross (
            ch_bam.target
                .filter { it -> it[0] != "" }
        )
        .map {
            row -> [ row[1][1][0], row[1][1][1], row[0][1][1] ]
        }
        .set { ch_bam_paired }
    //ch_bam_paired | view
    // EXAMPLE CHANNEL STRUCT: [[TARGET META], TARGET_BAM, CONTROL_BAM]


    /*
    * Call peaks with IgG
    */

    MACS2_NARROW_IGG (
        ch_bam_paired,
        params.macs_gsize
    )
    ch_narrow_igg       = MACS2_NARROW_IGG.out.peak
    ch_narrow_igg_summits     = MACS2_NARROW_IGG.out.bed
    ch_versions = MACS2_NARROW_IGG.out.versions

    MACS2_BROAD_IGG (
        ch_bam_paired,
        params.macs_gsize
    )
    ch_broad_igg       = MACS2_BROAD_IGG.out.peak
    ch_broad_igg_summits     = MACS2_BROAD_IGG.out.bed
    ch_versions = ch_versions.mix(MACS2_BROAD_IGG.out.versions)

    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //MACS2_BROAD_IGG.out.peak | view



    /*
    * Call peaks with no IgG
    */
    ch_bam.target.map{ it -> [ it[1][0], it[1][1], [] ] }
    .set { ch_bam_target_fctrl }
    //ch_bam_target_fctrl | view
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]

    MACS2_NARROW_NOIGG (
        ch_bam_target_fctrl,
        params.macs_gsize
    )

    MACS2_NARROW_NOIGG.out.peak
        .branch {
            with_control: it[0].control_group != ""
            no_control: it[0].control_group == ""
        }
        .set { ch_narrow_noigg }
    //ch_peaks_summits_noigg_narrow     = MACS2_NARROW_NOIGG.out.bed
    ch_versions = ch_versions.mix(MACS2_NARROW_NOIGG.out.versions)

    MACS2_BROAD_NOIGG (
        ch_bam_target_fctrl,
        params.macs_gsize
    )
    // MACS2_BROAD_NOIGG.out.peak | view
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    MACS2_BROAD_NOIGG.out.peak
        .branch {
            with_control: it[0].control_group != ""
            no_control: it[0].control_group == ""
        }
        .set { ch_broad_noigg }
    //ch_peaks_summits_noigg_broad     = MACS2_BROAD_NOIGG.out.bed
    ch_versions = ch_versions.mix(MACS2_BROAD_NOIGG.out.versions)


    /*
    * Pair igg and noigg MACS2 narrow peaks when available
    */
    NARROW_BEDTOOLS_INTERSECT(
            ch_narrow_igg
            .join ( ch_narrow_noigg.with_control.ifEmpty([[:],[]]) )
            .map { row -> [row[0], row[1], row[2]] },
            [[:],[]]
    )
    ch_narrow_filtered =  NARROW_BEDTOOLS_INTERSECT.out.intersect
    ch_versions = ch_versions.mix(NARROW_BEDTOOLS_INTERSECT.out.versions)
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.intersect | view

    NARROW_FILTER(
            ch_narrow_filtered
    )
    ch_narrow_filtered =  NARROW_FILTER.out.file
    ch_versions = ch_versions.mix(NARROW_FILTER.out.versions)

    /*
    * CHANNEL: pair igg and noigg MACS2 broad peaks
    */
    BROAD_BEDTOOLS_INTERSECT(
            ch_broad_noigg.with_control.ifEmpty([[:],[]])
            .join ( ch_broad_igg )
            .map {row -> [row[0], row[1], row[2]]},
            [[:],[]]
    )
    ch_broad_filtered =  BROAD_BEDTOOLS_INTERSECT.out.intersect
    ch_versions = ch_versions.mix(BROAD_BEDTOOLS_INTERSECT.out.versions)
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //BROAD_BEDTOOLS_INTERSECT.out.intersect | view

    BROAD_FILTER(
            ch_broad_filtered
    )
    ch_broad_filtered =  BROAD_FILTER.out.file
    ch_versions = ch_versions.mix(BROAD_FILTER.out.versions)


    /*
    * Call peaks against combined IgG
    */
    if (combine_igg){
        meta_combine_igg
            .map { it -> [ it.id, it ]}
            .set { ch_meta_comb_igg }

        MERGE_BAM_TARGET_MARKDUP (
            ch_meta_comb_igg
            .join(
                bam_markdup
                    .map { it -> [ it[0].id, it[1] ] } // [ val(meta.id), path(bam) ]
            )
            .filter { it -> it[1].is_control == true }
            .map { it -> [ it[1].control_group, it[2] ]} // [ val(meta.control_group), path(bam) ]
            .groupTuple()
        )

        MERGE_BAM_TARGET_MARKDUP.out.bam // [ val(control_group), path(control_bam) ]
            .cross(
                bam_markdup
                .filter { it -> it[0].is_control == false }
                .map    { it -> [ it[0].id, it[1] ] } // [ val(meta.id), path(target_bam) ]
                .join ( ch_meta_comb_igg ) // [ val(meta.id), path(target_bam), [meta_comb_igg] ]
                .map { it -> [ it[2].control_group, it[2], it[1] ]} // [ val(meta.control_group), [meta_comb_igg], path(target_bam) ]
            )
            .map { it -> [ it[1][1], it[1][2], it[0][1] ] } // EXAMPLE CHANNEL STRUCT: [[TARGET META], TARGET_BAM, CONTROL_BAM]
            .set { ch_bam_paired_comb_igg }



        MACS2_NARROW_COMB_IGG (
            ch_bam_paired_comb_igg,
            params.macs_gsize
        )
        ch_narrow_comb_igg       = MACS2_NARROW_COMB_IGG.out.peak
        ch_narrow_comb_igg_summits = MACS2_NARROW_COMB_IGG.out.bed
        ch_versions = ch_versions.mix(MACS2_NARROW_COMB_IGG.out.versions)

        NARROW_COMB_IGG_BEDTOOLS_INTERSECT(
            ch_narrow_comb_igg
            .join (
                ch_narrow_noigg.with_control
                    .map { it -> [ it[0].id, it[1] ]}
                    .join( ch_meta_comb_igg )
                    .map { it -> [ it[2], it[1] ] } // [ [meta], path(peak) ]
                .ifEmpty([[:],[]])
            ),
            [[:],[]]
        )
        ch_versions = ch_versions.mix(NARROW_COMB_IGG_BEDTOOLS_INTERSECT.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //NARROW_COMB_IGG_BEDTOOLS_INTERSECT.out.intersect | view

        NARROW_COMB_IGG_FILTER(
            NARROW_COMB_IGG_BEDTOOLS_INTERSECT.out.intersect
        )
        ch_narrow_comb_igg_filtered =  NARROW_COMB_IGG_FILTER.out.file
        ch_versions = ch_versions.mix(NARROW_COMB_IGG_FILTER.out.versions)


        MACS2_BROAD_COMB_IGG (
            ch_bam_paired_comb_igg,
            params.macs_gsize
        )
        ch_broad_comb_igg       = MACS2_BROAD_COMB_IGG.out.peak
        ch_broad_comb_igg_summits     = MACS2_BROAD_COMB_IGG.out.bed
        ch_versions = ch_versions.mix(MACS2_BROAD_COMB_IGG.out.versions)

        BROAD_COMB_IGG_BEDTOOLS_INTERSECT(
            ch_broad_comb_igg
            .join (
                ch_broad_noigg.with_control
                    .map { it -> [ it[0].id, it[1] ]}
                    .join( ch_meta_comb_igg )
                    .map { it -> [ it[2], it[1] ] } // [ [meta], path(peak) ]
                .ifEmpty([[:],[]])
            ),
            [[:],[]]
        )
        ch_versions = ch_versions.mix(BROAD_COMB_IGG_BEDTOOLS_INTERSECT.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //BROAD_COMB_IGG_BEDTOOLS_INTERSECT.out.intersect | view

        BROAD_COMB_IGG_FILTER(
            BROAD_COMB_IGG_BEDTOOLS_INTERSECT.out.intersect
        )
        ch_broad_comb_igg_filtered =  BROAD_COMB_IGG_FILTER.out.file
        ch_versions = ch_versions.mix(BROAD_COMB_IGG_FILTER.out.versions)

    }


    emit:
    versions = ch_versions
    narrow_igg = ch_narrow_igg
    narrow_noigg = MACS2_NARROW_NOIGG.out.peak
    narrow_filtered = ch_narrow_filtered
    narrow_comb_igg = ch_narrow_comb_igg
    narrow_comb_igg_filtered = ch_narrow_comb_igg_filtered

    broad_igg = ch_broad_igg
    broad_noigg = MACS2_BROAD_NOIGG.out.peak
    broad_filtered = ch_broad_filtered
    broad_comb_igg = ch_broad_comb_igg
    broad_comb_igg_filtered = ch_broad_comb_igg_filtered
}
