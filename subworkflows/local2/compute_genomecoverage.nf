/*
 * Convert bam files to bedgraph and bigwig with apropriate normalisation
 */

include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_MARKDUP_UNNORM   } from "../../modules/nf-core/bedtools/genomecov/main"
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_MARKDUP_SPIKEIN    } from "../../modules/nf-core/bedtools/genomecov/main"
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_DEDUP_UNNORM   } from "../../modules/nf-core/bedtools/genomecov/main"
//include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_DEDUP_SPIKEIN    } from "../../modules/nf-core/bedtools/genomecov/main"
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_MARKDUP } from "../../modules/local/for_patch/deeptools/bamcoverage/main"
include { DEEPTOOLS_BAMCOVERAGE as DEEPTOOLS_BAMCOVERAGE_DEDUP } from "../../modules/local/for_patch/deeptools/bamcoverage/main"

workflow COMPUTE_GENOMECOVERAGE {
    take:
    ch_bam_markdup         // channel: [ val(meta), [ bam ] ]
    ch_bai_markdup         // channel: [ val(meta), [ bai ] ]
    ch_bam_dedup         // channel: [ val(meta), [ bam ] ]
    ch_bai_dedup         // channel: [ val(meta), [ bai ] ]
    ch_chrom_sizes // channel: [ sizes ]
    ch_dummy_file  // channel: [ dummy ]
    norm_mode      // value:   ["Spikein", "RPKM", "CPM", "BPM", "RPGC", "None" ]
    metadata       // channel  [ csv ]

    main:
    ch_versions = Channel.empty()
    ch_bedgraph_markdup_unnorm = Channel.empty()
    ch_bedgraph_dedup_unnorm = Channel.empty()
    ch_bigwig_markdup = Channel.empty()
    ch_bigwig_dedup = Channel.empty()
    ch_bedgraph_spikein = Channel.empty()

    /*
    * Convert markdup.bam files to markdup.unnormalized.bedgraph
    * CHANNEL: Assign scale factor of 1
    */
    ch_bam_markdup
        .filter ( it -> it[0].is_control == false )
        .map { row -> [ row[0], row[1], 1 ] }
        .set { ch_bam_markdup_scale_factor_unnorm }
        //ch_bam_markdup_scale_factor_unnorm | view

    BEDTOOLS_GENOMECOV_MARKDUP_UNNORM (
        ch_bam_markdup_scale_factor_unnorm,
        ch_dummy_file,
        "bedGraph"
    )

    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_MARKDUP_UNNORM.out.versions)
    ch_bedgraph_markdup_unnorm = BEDTOOLS_GENOMECOV_MARKDUP_UNNORM.out.genomecov
    //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
    //BEDTOOLS_GENOMECOV_MARKDUP_UNNORM.out.genomecov | view

    /*
    * Convert dedup.bam files to dedup.unnormalized.bedgraph
    * CHANNEL: Assign scale factor of 1
    */
    ch_bam_dedup
        .filter ( it -> it[0].is_control == true )
        .map { row -> [ row[0], row[1], 1 ] }
        .set { ch_bam_dedup_scale_factor_unnorm }
        //ch_bam_dedup_scale_factor_unnorm | view

    BEDTOOLS_GENOMECOV_DEDUP_UNNORM (
        ch_bam_dedup_scale_factor_unnorm,
        ch_dummy_file,
        "bedGraph"
    )

    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_DEDUP_UNNORM.out.versions)
    ch_bedgraph_dedup_unnorm = BEDTOOLS_GENOMECOV_DEDUP_UNNORM.out.genomecov
    //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
    //BEDTOOLS_GENOMECOV_DEDUP_UNNORM.out.genomecov | view

    /*
    * Convert markdup.bam files to markdup.CPM.bigwig for genome browser
    */
    /*
     * CHANNEL: Combine bam and bai files on id
    */
    ch_bam_markdup
        .map { row -> [row[0].id, row ].flatten()}
        .join ( ch_bai_markdup.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4], 1] }
        .set { ch_bam_bai_markdup_scale_factor }
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
    //ch_bam_bai_markdup_scale_factor | view

    /*
    * MODULE: Convert bam files to bedgraph and normalise
    */
    DEEPTOOLS_BAMCOVERAGE_MARKDUP (
        ch_bam_bai_markdup_scale_factor
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_MARKDUP.out.versions)
    ch_bigwig_markdup = DEEPTOOLS_BAMCOVERAGE_MARKDUP.out.bigwig
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
    //ch_bedgraph | view


    /*
    * Convert dedup.bam files to dedup.CPM.bigwig for genome browser
    */
    /*
     * CHANNEL: Combine bam and bai files on id
    */
    ch_bam_dedup
        .map { row -> [row[0].id, row ].flatten()}
        .join ( ch_bai_dedup.map { row -> [row[0].id, row ].flatten()} )
        .map { row -> [row[1], row[2], row[4], 1] }
        .filter ( it -> it[0].is_control == true )
        .set { ch_bam_bai_dedup_scale_factor }
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI, SCALE_FACTOR]
    //ch_bam_bai_dedup_scale_factor | view

    /*
    * MODULE: Convert bam files to bedgraph and normalise
    */
    DEEPTOOLS_BAMCOVERAGE_DEDUP (
        ch_bam_bai_dedup_scale_factor
    )
    ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE_DEDUP.out.versions)
    ch_bigwig_dedup = DEEPTOOLS_BAMCOVERAGE_DEDUP.out.bigwig
    // EXAMPLE CHANNEL STRUCT: [[META], BAM, BAI]
    //ch_bedgraph_dedup | view


    emit:
    bedgraph_markdup_unnorm = ch_bedgraph_markdup_unnorm        // channel: [ val(meta), [ bedgraph ] ]
    bedgraph_dedup_unnorm = ch_bedgraph_dedup_unnorm        // channel: [ val(meta), [ bedgraph ] ]
    bigwig_markdup = ch_bigwig_markdup // channel: [val(meta), bigwig]
    bigwig_dedup = ch_bigwig_dedup // channel: [val(meta), bigwig]
    versions = ch_versions                      // channel: [ versions.yml ]
}
