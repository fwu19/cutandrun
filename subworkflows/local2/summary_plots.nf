/*
 * Convert bam files to bedgraph and bigwig with apropriate normalisation
 */

include { GTF2GENES                                 } from "../../modules/local2/gtf2genes"
include { BIGWIG_AVERAGE                            } from "../../modules/local2/bigwig_average"
include { BIGWIG_AVERAGE as BIGWIG_AVERAGE_IGG      } from "../../modules/local2/bigwig_average"
include { TORNADO_PLOTS                             } from "../../modules/local2/tornado_plots"
include { UPDATE_PEAKS                              } from "../../modules/local2/update_peaks"
include { TORNADO_PLOTS  as TORNADO_PLOTS_COMB_IGG              } from "../../modules/local2/tornado_plots"
include { UPDATE_PEAKS as UPDATE_PEAKS_COMB_IGG                 } from "../../modules/local2/update_peaks"

workflow SUMMARY_PLOTS {
    take:
    ch_bigwig         // channel: [ val(meta), [ bam ] ]
    ch_conp_bed       // [ target, path(conp) ]
    ch_conp_ann       // [ target, path(conp) ]
    ch_dp
    ch_conp_bed_comb_igg
    ch_conp_ann_comb_igg
    ch_dp_comb_igg
    use_igg
    gene_gtf
    gene_bed

    main:
    ch_versions = Channel.empty()
    ch_avg_bigwig = Channel.empty()
    ch_versions = Channel.empty()

    /*
    * convert gtf to gene.bed
    */
    if (gene_bed =~ 'dummy'){
        GTF2GENES(gene_gtf)
        gene_bed = GTF2GENES.out.bed
    }

    /*
    * Convert markdup.bam files to markdup.CPM.bigwig for genome browser
    */

    BIGWIG_AVERAGE(
        ch_bigwig
            .filter { it[0].is_control == false }
            .map { it -> [ [ it[0].sample_group, it[0].target ], it[1] ] }
            .groupTuple()
    )
    ch_versions = ch_versions.mix(BIGWIG_AVERAGE.out.versions)


    ch_avg_bigwig = BIGWIG_AVERAGE.out.bigwig
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
    //ch_bedgraph_dedup | view

    if ('individual' in use_igg){
        BIGWIG_AVERAGE_IGG(
            ch_bigwig
            .filter { it[0].is_control == true }
            .map { it -> [ [ it[0].sample_group, it[0].target ], it[1] ] }
            .groupTuple()
        // [[sample_group, target], Bigwig1, Bigwig2, ...]
        )
        ch_versions = ch_versions.mix(BIGWIG_AVERAGE_IGG.out.versions)

        /*
        * Update peak bed files by adding peak annotation, differential peaks and so on
        */
        UPDATE_PEAKS(
            ch_conp_bed.collect{it[1]},
            ch_conp_ann.collect(),
            ch_dp
        )
        ch_versions = ch_versions.mix(UPDATE_PEAKS.out.versions)

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
            .map { it -> [ it[1][0], it[1][1], it[0][1] ] }, // [ val(target), [bigwig], [conp_bed] ]
            gene_bed
        )
        ch_versions = ch_versions.mix(TORNADO_PLOTS.out.versions)

    }

    if ([ 'group', 'all', 'custom' ].any { it in use_igg }){
        /*
        * Update peak bed files by adding peak annotation, differential peaks and so on
        */
        UPDATE_PEAKS_COMB_IGG(
            ch_conp_bed_comb_igg.collect{it[1]},
            ch_conp_ann_comb_igg.collect(),
            ch_dp_comb_igg
        )
        ch_versions = ch_versions.mix(UPDATE_PEAKS_COMB_IGG.out.versions)

        /*
        * tornado plots for each antibody
        */
        TORNADO_PLOTS_COMB_IGG(
            ch_conp_bed_comb_igg
            .cross (
                ch_avg_bigwig
                    .map { it -> [ it[0][1], it[1] ] }
                    .groupTuple()
            )
            .map { it -> [ it[1][0], it[1][1], it[0][1] ] }, // [ val(target), [bigwig], [conp_bed] ]
            gene_bed
        )
        ch_versions = ch_versions.mix(TORNADO_PLOTS_COMB_IGG.out.versions)

    }

    emit:
    versions = ch_versions                      // channel: [ versions.yml ]
}
