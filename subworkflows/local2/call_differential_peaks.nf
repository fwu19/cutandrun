/*
* Call peaks with both MACS2 and SEACR
*/

include { READS_IN_CONSENSUS_PEAKS                             } from '../../modules/local2/reads_in_consensus_peaks'
include { COLLECT_COUNT_MATRIX                                 } from '../../modules/local2/collect_count_matrix'
include { DIFFERENTIAL_PEAKS                                   } from '../../modules/local2/differential_peaks'
include { READS_IN_CONSENSUS_PEAKS as READS_IN_CONSENSUS_PEAKS_COMB_IGG     } from '../../modules/local2/reads_in_consensus_peaks'
include { COLLECT_COUNT_MATRIX as COLLECT_COUNT_MATRIX_COMB_IGG             } from '../../modules/local2/collect_count_matrix'
include { DIFFERENTIAL_PEAKS as DIFFERENTIAL_PEAKS_COMB_IGG                 } from '../../modules/local2/differential_peaks'


workflow CALL_DIFFERENTIAL_PEAKS {
    take:
    samplesheet
    conp_bed
    samplesheet_combine_igg
    conp_bed_comb_igg
    samtools_bam
    use_igg
    srcdir


    main:
    ch_cts = Channel.empty()
    ch_conp_counts = Channel.empty()
    ch_dp = Channel.empty()
    ch_cts_comb_igg = Channel.empty()
    ch_conp_counts_comb_igg = Channel.empty()
    ch_dp_comb_igg = Channel.empty()
    ch_versions = Channel.empty()

    // get pre-generated alignment files
    if (!params.run_alignment){
        samtools_bam = Channel.fromPath("${srcdir}/csv/map2genome.${params.aligner}.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it , file("${srcdir}/${it.target_bam}", checkIfExists: true) ] }
    }


    if ('individual' in use_igg){
        // get pre-generated peak files
        if ( !params.run_peak_qc ) {
            conp_bed = Channel.fromPath("${srcdir}/csv/consensus_peaks.individual_igg.csv", checkIfExists: true)
                    .splitCsv(header: true)
                    .map { it -> [ it.target, file("${srcdir}/${it.conp_bed}", checkIfExists: true) ]}
                    .groupTuple(by: 0)
        }

        // Count reads in consensus peaks
        READS_IN_CONSENSUS_PEAKS(
            conp_bed
                .cross (
                    samtools_bam
                            .filter { it[0].target != "IgG" }
                            .map { meta -> [ meta[0].target, meta ] }
                )
                .map { bed, bam -> [ bam[1][0], bam[1][1], bed[1] ] }
        )
        ch_versions = ch_versions.mix(READS_IN_CONSENSUS_PEAKS.out.versions)
        ch_conp_counts = READS_IN_CONSENSUS_PEAKS.out.count

        COLLECT_COUNT_MATRIX(
            ch_conp_counts
                .groupTuple( by: 0 )
                .map { meta, counts -> tuple(meta, counts.flatten())}
                .cross ( conp_bed )
                .map { counts, bed -> tuple(counts[0], counts[1], bed[1]) }
                .combine(samplesheet)
        )
        ch_cts = COLLECT_COUNT_MATRIX.out.cts
        ch_versions = ch_versions.mix(COLLECT_COUNT_MATRIX.out.versions)
        }

        DIFFERENTIAL_PEAKS(
                ch_cts
                .combine(samplesheet)
                .combine(Channel.fromPath( params.comparison, checkIfExists: true ))
        )
        ch_dp = DIFFERENTIAL_PEAKS.out.data
        ch_versions = ch_versions.mix(DIFFERENTIAL_PEAKS.out.versions)


    if ([ 'group', 'all', 'custom' ].any { it in use_igg }){
        // get pre-generated peak files
        if ( !params.run_peak_qc) {
            conp_bed_comb_igg = Channel.fromPath("${srcdir}/csv/consensus_peaks.combined_igg.csv")
                    .splitCsv(header: true)
                    .map { it -> [ it.target, file("${srcdir}/${it.conp_bed}", checkIfExists: true) ]}
                    .groupTuple(by: 0)
        }

        // Count reads in consensus peaks
        READS_IN_CONSENSUS_PEAKS_COMB_IGG(
            conp_bed_comb_igg
                .cross (
                    samtools_bam
                            .filter { it[0].target != "IgG" }
                            .map { meta -> [ meta[0].target, meta ] }
                )
                .map { bed, bam -> [ bam[1][0], bam[1][1], bed[1] ] }
        )
        ch_versions = ch_versions.mix(READS_IN_CONSENSUS_PEAKS_COMB_IGG.out.versions)
        ch_conp_counts_comb_igg = READS_IN_CONSENSUS_PEAKS_COMB_IGG.out.count

        COLLECT_COUNT_MATRIX_COMB_IGG(
            ch_conp_counts_comb_igg
                .groupTuple( by: 0 )
                .map { meta, counts -> tuple(meta, counts.flatten())}
                .cross ( conp_bed_comb_igg )
                .map { counts, bed -> tuple(counts[0], counts[1], bed[1]) }
                .combine(samplesheet_combine_igg)
        )
        ch_cts_comb_igg = COLLECT_COUNT_MATRIX_COMB_IGG.out.cts
        ch_versions = ch_versions.mix(COLLECT_COUNT_MATRIX_COMB_IGG.out.versions)

        DIFFERENTIAL_PEAKS_COMB_IGG(
                ch_cts_comb_igg
                .combine(samplesheet_combine_igg)
                .combine(Channel.fromPath( params.comparison, checkIfExists: true ))
        )
        ch_dp_comb_igg = DIFFERENTIAL_PEAKS_COMB_IGG.out.data
        ch_versions = ch_versions.mix(DIFFERENTIAL_PEAKS_COMB_IGG.out.versions)
    }



    emit:
    versions = ch_versions
    dp = ch_dp
    dp_comb_igg = ch_dp_comb_igg

}
