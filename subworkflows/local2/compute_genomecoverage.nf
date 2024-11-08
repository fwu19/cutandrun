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
    ch_bam_markdup.map { row ->
            [ row[0], row[1], 1 ]
        }
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
    ch_bam_dedup.map { row ->
            [ row[0], row[1], 1 ]
        }
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
    * Convert markdup.bam files to spike-in normalized bedgraph if required
    */

    if (norm_mode == "Spikein") {
        /*
        * CHANNEL: Load up alignment metadata into channel
        */
        metadata.splitCsv ( header:true, sep:"," )
            .map { row -> [ row[0].id, row[1] ]}
            .set { ch_metadata }
        //ch_metadata | view

        /*
        * CHANNEL: Calculate scale factor for each sample based on a constant devided by the number
        *          of reads aligned to the spike-in genome.
        */
        ch_bam_markdup.map { row -> [ row[0].id, row[0], row[1] ]}
            .join ( ch_metadata )
            .map { row ->
                def denominator = row[3].find{ it.key == "bt2_total_aligned" }?.value.toInteger()
                [ row[1], row[2], params.normalisation_c / (denominator != 0 ? denominator : params.normalisation_c) ]
            }
            .set { ch_bam_markdup_scale_factor_spikein }
        // EXAMPLE CHANNEL STRUCT: [id, scale_factor]
        //ch_bam_markdup_scale_factor | view


        /*
        * MODULE: Convert bam files to bedgraph
        */
        BEDTOOLS_GENOMECOV_MARKDUP_SPIKEIN (
            ch_bam_markdup_scale_factor_spikein,
            ch_dummy_file,
            "markdup.spikein_norm.bedGraph"
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV_MARKDUP_SPIKEIN.out.versions)
        ch_bedgraph_spikein = BEDTOOLS_GENOMECOV_MARKDUP_SPIKEIN.out.genomecov
        //EXAMPLE CHANNEL STRUCT: [META], BEDGRAPH]
        //BEDTOOLS_GENOMECOV_MARKDUP_SPIKEIN.out.genomecov | view

        /*
        * CHANNEL: Dump scale factor values
        */
        if(params.dump_scale_factors) {
            ch_scale_factor = ch_bam_markdup_scale_factor_spikein
            .map { [it[0].id, it[0].group, it[2]] }
            .toSortedList( { a, b -> a[0] <=> b[0] } )
            .map { list ->
                new File('scale-factors.markdup.by-spikein.csv').withWriter('UTF-8') { writer ->
                    list.each { item ->
                        str = item[0] + "," + item[1] + "," + item[2]
                        writer.write(str + "\n")
                    }
                }
            }
        }

    }


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
    bedgraph_spikein   = ch_bedgraph_spikein // channel: [ val(meta), [ bedgraph ] ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
