/*
 * Convert bam files to bedgraph and bigwig with apropriate normalisation
 */

include { BIGWIG_AVERAGE                            } from "../../modules/local2/bigwig_average"
include { BIGWIG_AVERAGE as BIGWIG_AVERAGE_IGG } from "../../modules/local2/bigwig_average"
include { TORNADO_PLOTS                             } from "../../modules/local2/tornado_plots"
include { UPDATE_PEAKS                              } from "../../modules/local2/update_peaks"

workflow SUMMARY_PLOTS {
    take:
    ch_bigwig         // channel: [ val(meta), [ bam ] ]
    ch_conp_bed
    ch_conp_ann
    ch_dp
    gene_bed
    average_igg

    main:
    ch_versions = Channel.empty()
    ch_avg_bigwig = Channel.empty()


    /*
    * Convert markdup.bam files to markdup.CPM.bigwig for genome browser
    */

    BIGWIG_AVERAGE(
        ch_bigwig
            .filter { it[0].is_control == false }
            .map { it -> [ [ it[0].sample_group, it[0].target ], it[1] ] }
            .groupTuple()
    )

    if (average_igg){
        BIGWIG_AVERAGE_IGG(
            ch_bigwig
            .filter { it[0].is_control == true }
            .map { it -> [ [ it[0].sample_group, it[0].target ], it[1] ] }
            .groupTuple()
        // [[sample_group, target], Bigwig1, Bigwig2, ...]
        )
    }

    ch_versions = BIGWIG_AVERAGE.out.versions
    ch_avg_bigwig = BIGWIG_AVERAGE.out.bigwig
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
    //ch_bedgraph_dedup | view

    /*
    * Update peak bed files by adding peak annotation, differential peaks and so on
    */
    UPDATE_PEAKS(
        ch_conp_bed.collect{it[1]},
        ch_conp_ann,
        ch_dp
    )

    /*
    * tornado plots for each antibody
    */
    TORNADO_PLOTS(
        ch_conp_bed
            .cross (
                ch_avg_bigwig
                    .map { it -> [ it[0][1], it[1] ] }
                    .groupTuple()
            )
            .map { it -> [ it[1][0], it[1][1], it[0][1] ] },
        gene_bed
        // [ val(target), [bigwig], [conp_bed] ]
    )


    emit:
    bigwig = ch_avg_bigwig        // channel: [ val(group), [ bedgraph ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
