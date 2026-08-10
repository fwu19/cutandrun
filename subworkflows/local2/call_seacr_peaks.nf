/*
* Call peaks with SEACR
*/

include { SEACR_CALLPEAK as SEACR_CALLPEAK_IGG                          } from "../../modules/nf-core/seacr/callpeak/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_NOIGG                        } from "../../modules/nf-core/seacr/callpeak/main"
include { BEDTOOLS_INTERSECT as SEACR_BEDTOOLS_INTERSECT                } from "../../modules/nf-core/bedtools/intersect/main"
include { UNION_BEDGRAPH as UNION_BEDGRAPH_DEDUP                        } from "../../modules/local2/union_bedgraph"
include { MERGE_UNION_BEDGRAPH as MERGE_UNION_BEDGRAPH_DEDUP            } from "../../modules/local2/merge_union_bedgraph"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_COMB_IGG                     } from "../../modules/nf-core/seacr/callpeak/main"
include { BEDTOOLS_INTERSECT as SEACR_BEDTOOLS_INTERSECT_COMB_IGG       } from "../../modules/nf-core/bedtools/intersect/main"

workflow CALL_SEACR_PEAKS {
    take:
    bedgraph_markdup
    bedgraph_dedup
    combine_igg
    meta_combine_igg
    skip_individual_igg
    srcdir

    main:

    ch_seacr_igg                = Channel.empty()
    ch_seacr_noigg              = Channel.empty()
    ch_seacr_filtered           = Channel.empty()
    ch_seacr_comb_igg           = Channel.empty()
    ch_seacr_comb_igg_filtered  = Channel.empty()

    if (!params.run_alignment){
        Channel.fromPath("${srcdir}/csv/genome_coverage.markdup_bdg.csv")
            .splitCsv(header: true)
            .map {row ->
                row.is_control = row.is_control.toBoolean()
                def meta_map = row.clone()
                def column_names = meta_map.keySet().toList()
                def last_column_name = column_names[-1]
                def last_column_value = meta_map.remove(last_column_name)
                def prefixed_value = "${srcdir}/${last_column_value}"
                return [ meta_map, prefixed_value ]
            }
            .set{ bedgraph_markdup }

        Channel.fromPath("${srcdir}/csv/genome_coverage.dedup_bdg.csv")
            .splitCsv(header: true)
            .map {row ->
                row.is_control = row.is_control.toBoolean()
                def meta_map = row.clone()
                def column_names = meta_map.keySet().toList()
                def last_column_name = column_names[-1]
                def last_column_value = meta_map.remove(last_column_name)
                def prefixed_value = "${srcdir}/${last_column_value}"
                return [ meta_map, prefixed_value ]
            }
            .set {bedgraph_dedup }
    }

    bedgraph_markdup
        .filter { it -> it[0].is_control == false }
        .map{
            it -> [ it[0].control_group, it[0], it[1] ]
        }
        .set { ch_bedgraph_target }
    //ch_bedgraph_target | view

    bedgraph_dedup
        .filter { it -> it[0].is_control == true }
        .map{
            it -> [ it[0].control_group, it[0], it[1] ]
        }
        .set{ ch_bedgraph_control }
    //ch_bedgraph_control | view

    ch_bedgraph_control
        .cross(
            ch_bedgraph_target
                .filter { it -> it[0] != "" }
        )
        .map { it -> [ it[1][1], it[1][2], it[0][2] ] }
        .set { ch_bedgraph_paired }
    //println "bedgraph_paired"
    //ch_bedgraph_paired | view
    // EXAMPLE CHANNEL STRUCT: [[TARGET META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]

    /*
    * without IgG control
    */
    ch_bedgraph_target
        .map{ it -> [ it[1], it[2], [] ] }
        .set { ch_bedgraph_target_fctrl }
    // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
    //println "bedgraph_target_fctrl"
    //ch_bedgraph_target_fctrl | view

    SEACR_CALLPEAK_NOIGG (
        ch_bedgraph_target_fctrl,
        params.seacr_peak_threshold
    )
    SEACR_CALLPEAK_NOIGG.out.bed
        .branch {
            with_control: it[0].control_group != ""
            no_control: it[0].control_group == ""
        }
        .set { ch_seacr_noigg}
    ch_versions = SEACR_CALLPEAK_NOIGG.out.versions
    // EXAMPLE CHANNEL STRUCT: [[META], BED]
    //SEACR_NO_IGG.out.bed | view


    /*
    * with IgG control
    */

    if (!skip_individual_igg){
        SEACR_CALLPEAK_IGG (
            ch_bedgraph_paired,
            params.seacr_peak_threshold
        )
        ch_seacr_igg    = SEACR_CALLPEAK_IGG.out.bed
        ch_versions = ch_versions.mix(SEACR_CALLPEAK_IGG.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //SEACR_CALLPEAK_IGG.out.bed | view


        // intersect igg and no-igg peaks
        SEACR_BEDTOOLS_INTERSECT(
            ch_seacr_igg
            .join ( ch_seacr_noigg.with_control.ifEmpty([[:],[]]) )
            .map { it -> [it[0], it[1], it[2]]},
            [[:],[]]
        )
        ch_seacr_filtered = SEACR_BEDTOOLS_INTERSECT.out.intersect
        ch_versions = ch_versions.mix(SEACR_BEDTOOLS_INTERSECT.out.versions)
        // EXAMPLE CHANNEL STRUCT: [[META], BED]
        //BEDTOOLS_INTERSECT.out.intersect | view

    }

    /*
    * Call peaks against combined IgG
    */
    if (combine_igg){
        meta_combine_igg
            .map { it -> [ it.id, it ]}
            .set { ch_meta_comb_igg }

        UNION_BEDGRAPH_DEDUP (
            ch_meta_comb_igg
            .join(
                bedgraph_dedup
                    .map { it -> [ it[0].id, it[1] ] } // [ val(meta.id), path(bedgraph) ]
            )
            .filter { it -> it[1].is_control == true }
            .map { it -> [ it[1].control_group, it[2] ]} // [ val(meta.control_group), path(bedgraph) ]
            .groupTuple()
        )

        MERGE_UNION_BEDGRAPH_DEDUP (
            UNION_BEDGRAPH_DEDUP.out.bedgraph
        )

        MERGE_UNION_BEDGRAPH_DEDUP.out.bedgraph // [ val(control_group), path(control_bedgraph) ]
            .cross(
                bedgraph_markdup
                .filter { it -> it[0].is_control == false }
                .map    { it -> [ it[0].id, it[1] ] } // [ val(meta.id), path(target_bedgraph) ]
                .join ( ch_meta_comb_igg ) // [ val(meta.id), path(target_bedgraph), [meta_comb_igg] ]
                .map { it -> [ it[2].control_group, it[2], it[1] ]} // [ val(meta.control_group), [meta_comb_igg], path(target_bedgraph) ]
            )
            .map { it -> [ it[1][1], it[1][2], it[0][1] ] } // EXAMPLE CHANNEL STRUCT: [[TARGET META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]
            .set { ch_bedgraph_paired_comb_igg }


        SEACR_CALLPEAK_COMB_IGG (
            ch_bedgraph_paired_comb_igg,
            params.seacr_peak_threshold
        )
        ch_seacr_comb_igg = SEACR_CALLPEAK_COMB_IGG.out.bed
        ch_versions = SEACR_CALLPEAK_COMB_IGG.out.versions

        SEACR_BEDTOOLS_INTERSECT_COMB_IGG (
            ch_seacr_comb_igg // [ [meta], path(bed) ]
            .join (
                ch_seacr_noigg.with_control
                    .map { it -> [ it[0].id, it[1] ]}
                    .join( ch_meta_comb_igg )
                    .map { it -> [ it[2], it[1] ] } // [ [meta], path(bed) ]
                    .ifEmpty([[:],[]])
            ),
            [[:],[]]
        )
        ch_seacr_comb_igg_filtered = SEACR_BEDTOOLS_INTERSECT_COMB_IGG.out.intersect
        ch_versions = ch_versions.mix(SEACR_BEDTOOLS_INTERSECT_COMB_IGG.out.versions)

    }



    emit:
    versions = ch_versions
    seacr_noigg = SEACR_CALLPEAK_NOIGG.out.bed
    seacr_igg = ch_seacr_igg
    seacr_filtered = ch_seacr_filtered
    seacr_comb_igg = ch_seacr_comb_igg
    seacr_comb_igg_filtered = ch_seacr_comb_igg_filtered
}
