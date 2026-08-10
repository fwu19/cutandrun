nextflow.enable.dsl=2

/*
* Writes out csv files containing output paths
*/

include { WRITE_CSV as WRITE_CSV_BT2 } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MARKDUP_BDG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_DEDUP_BDG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_NARROW_NOIGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_NARROW_FILTERED } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_NARROW_IGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_NARROW_COMB_IGG_FILTERED } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_NARROW_COMB_IGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_BROAD_NOIGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_BROAD_FILTERED } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_BROAD_IGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_BROAD_COMB_IGG_FILTERED } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_MACS2_BROAD_COMB_IGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SEACR_NOIGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SEACR_FILTERED } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SEACR_IGG } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SEACR_COMB_IGG_FILTERED } from '../../modules/local2/write_csv.nf'
include { WRITE_CSV as WRITE_CSV_SEACR_COMB_IGG } from '../../modules/local2/write_csv.nf'


workflow WRITE_OUTPUT_CSV {
    take:
    samtools_bam_bai_markdup
    bedgraph_markdup
    bedgraph_dedup
    macs2_narrow_noigg
    macs2_narrow_filtered
    macs2_narrow
    macs2_narrow_comb_igg_filtered
    macs2_narrow_comb_igg
    macs2_broad_noigg
    macs2_broad_filtered
    macs2_broad
    macs2_broad_comb_igg_filtered
    macs2_broad_comb_igg
    seacr_noigg
    seacr_filtered
    seacr_igg
    seacr_comb_igg_filtered
    seacr_comb_igg

    main:

    /* map2genome.csv */
    if (params.run_alignment && params.workflow == "cutandrun" ){
        samtools_bam_bai_markdup
                .map {
                    it -> it[0] + [target_bam: "1_individual_samples/02_alignment/bowtie2/${params.target_genome}/${it[1].name}" ]
                    }
                .collect()
                .set{bt2}
        WRITE_CSV_BT2(
            bt2,
            "map2genome.bowtie2.csv"
        )

    }

    /* genome_coverage_bdg.csv */
    if(params.run_mark_dups && params.run_remove_dups){
        bedgraph_markdup
                .map {
                    it -> it[0] + [target_bdg: "1_individual_samples/03_genome_coverage/unnormalized_bedgraph/${it[1]?.name ?: ' '}"]
                    }
                .collect()
                .set{bdg_md}
        WRITE_CSV_MARKDUP_BDG(
            bdg_md,
            "genome_coverage.markdup_bdg.csv"
        )

        bedgraph_dedup
                .map {
                    it -> it[0] + [target_bdg: "1_individual_samples/03_genome_coverage/unnormalized_bedgraph/${it[1]?.name ?: ' '}"]
                    }
                .collect()
                .set{bdg_dd}
        WRITE_CSV_DEDUP_BDG(
            bdg_dd,
            "genome_coverage.dedup_bdg.csv"
        )

    }

    /* peaks.csv */
    if (params.run_peak_calling && params.workflow == "cutandrun" ){
        macs2_narrow_noigg
                .map {
                    it -> it[0] + [macs2_narrow_noigg: "${params.macs2_narrow_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_narrow_noigg}
        WRITE_CSV_MACS2_NARROW_NOIGG(
            macs2_narrow_noigg,
                "peaks.macs2_narrow_noigg.csv"
        )

        macs2_narrow_filtered
                .map {
                    it -> it[0] + [macs2_narrow_filtered: "${params.macs2_narrow_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_narrow_filtered}
        WRITE_CSV_MACS2_NARROW_FILTERED(
            macs2_narrow_filtered,
                "peaks.individual_igg.macs2_narrow_filtered.csv"
        )

        macs2_narrow
                .map {
                    it -> it[0] + [macs2_narrow_igg: "${params.macs2_narrow_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_narrow_igg}
        WRITE_CSV_MACS2_NARROW_IGG(
            macs2_narrow_igg,
                "peaks.individual_igg.macs2_narrow_igg.csv"
        )

        macs2_narrow_comb_igg_filtered
                .map {
                    it -> it[0] + [macs2_narrow_filtered: "${params.macs2_narrow_comb_igg_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_narrow_comb_igg_filtered}
        WRITE_CSV_MACS2_NARROW_COMB_IGG_FILTERED(
            macs2_narrow_comb_igg_filtered,
                "peaks.combine_igg.macs2_narrow_filtered.csv"
        )

        macs2_narrow_comb_igg
                .map {
                    it -> it[0] + [macs2_narrow_igg: "${params.macs2_narrow_comb_igg_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_narrow_comb_igg}
        WRITE_CSV_MACS2_NARROW_COMB_IGG(
            macs2_narrow_comb_igg,
                "peaks.combine_igg.macs2_narrow_igg.csv"
        )

    }

    if (params.run_peak_calling && params.workflow == "cutandrun" ){
        macs2_broad_noigg
                .map {
                    it -> it[0] + [macs2_broad_noigg: "${params.macs2_broad_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_broad_noigg}
        WRITE_CSV_MACS2_BROAD_NOIGG(
            macs2_broad_noigg,
                "peaks.macs2_broad_noigg.csv"
        )

        macs2_broad_filtered
                .map {
                    it -> it[0] + [macs2_broad_filtered: "${params.macs2_broad_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_broad_filtered}
        WRITE_CSV_MACS2_BROAD_FILTERED(
            macs2_broad_filtered,
                "peaks.individual_igg.macs2_broad_filtered.csv"
        )

        macs2_broad
                .map {
                    it -> it[0] + [macs2_broad_igg: "${params.macs2_broad_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_broad_igg}
        WRITE_CSV_MACS2_BROAD_IGG(
            macs2_broad_igg,
                "peaks.individual_igg.macs2_broad_igg.csv"
        )

        macs2_broad_comb_igg_filtered
                .map {
                    it -> it[0] + [macs2_broad_filtered: "${params.macs2_broad_comb_igg_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_broad_comb_igg_filtered}
        WRITE_CSV_MACS2_BROAD_COMB_IGG_FILTERED(
            macs2_broad_comb_igg_filtered,
                "peaks.combine_igg.macs2_broad_filtered.csv"
        )

        macs2_broad_comb_igg
                .map {
                    it -> it[0] + [macs2_broad_igg: "${params.macs2_broad_comb_igg_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{macs2_broad_comb_igg}
        WRITE_CSV_MACS2_BROAD_COMB_IGG(
            macs2_broad_comb_igg,
                "peaks.combine_igg.macs2_broad_igg.csv"
        )

    }

    if (params.run_peak_calling && params.workflow == "cutandrun" ){
        seacr_noigg
                .map {
                    it -> it[0] + [seacr_noigg: "${params.seacr_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{seacr_noigg}
        WRITE_CSV_SEACR_NOIGG(
            seacr_noigg,
                "peaks.seacr_noigg.csv"
        )

        seacr_filtered
                .map {
                    it -> it[0] + [seacr_filtered: "${params.seacr_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{seacr_filtered}
        WRITE_CSV_SEACR_FILTERED(
            seacr_filtered,
                "peaks.individual_igg.seacr_filtered.csv"
        )

        seacr_igg
                .map {
                    it -> it[0] + [seacr_igg: "${params.seacr_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{seacr_igg}
        WRITE_CSV_SEACR_IGG(
            seacr_igg,
                "peaks.individual_igg.seacr_igg.csv"
        )

        seacr_comb_igg_filtered
                .map {
                    it -> it[0] + [seacr_filtered: "${params.seacr_comb_igg_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{seacr_comb_igg_filtered}
        WRITE_CSV_SEACR_COMB_IGG_FILTERED(
            seacr_comb_igg_filtered,
                "peaks.combine_igg.seacr_filtered.csv"
        )

        seacr_comb_igg
                .map {
                    it -> it[0] + [seacr_igg: "${params.seacr_comb_igg_dir}/${it[1]?.name ?: ' '}"]
                }
                .collect()
                .set{seacr_comb_igg}
        WRITE_CSV_SEACR_COMB_IGG(
            seacr_comb_igg,
                "peaks.combine_igg.seacr_igg.csv"
        )

    }

}
