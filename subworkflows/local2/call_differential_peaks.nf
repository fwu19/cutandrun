/*
* Call peaks with both MACS2 and SEACR
*/

include { READS_IN_CONSENSUS_PEAKS                             } from '../../modules/local2/reads_in_consensus_peaks'
include { COLLECT_COUNT_MATRIX                                 } from '../../modules/local2/collect_count_matrix'
include { DIFFERENTIAL_PEAKS                                   } from '../../modules/local2/differential_peaks'


workflow CALL_DIFFERENTIAL_PEAKS {
    take:
    samplesheet
    conp_bed
    samtools_bam
    srcdir


    main:
        ch_cts = Channel.empty()
        ch_conp_counts = Channel.empty()
        ch_versions = Channel.empty()

        // Count reads in consensus peaks
        READS_IN_CONSENSUS_PEAKS(
            conp_bed
                .cross (
                    samtools_bam
                            .filter { it[0].target != "IgG" }
                            .map { it -> [ it[0].target, it ] }
                )
                .map { it -> [ it[1][1][0], it[1][1][1], it[0][1] ] }
        )
        ch_versions = ch_versions.mix(READS_IN_CONSENSUS_PEAKS.out.versions)
        ch_conp_counts = READS_IN_CONSENSUS_PEAKS.out.count

        COLLECT_COUNT_MATRIX(
            ch_conp_counts
                .groupTuple( by: 0 )
                .map { it -> [ it[0], it[1].flatten().collect() ]}
                .cross ( conp_bed )
                .map { it -> [ it[0][0], it[0][1], it[1][1] ]}
                .combine(samplesheet)
        )
        ch_cts = COLLECT_COUNT_MATRIX.out.cts
        ch_versions = ch_versions.mix(COLLECT_COUNT_MATRIX.out.versions)

        // if --comparison is a dummy file or empty file is used, write an error message and continue.
        // if file has contents but not a correct format, throw an error.
        if ( params.run_differential_peaks ){
            DIFFERENTIAL_PEAKS(
                ch_cts
                .combine(samplesheet)
                .combine(Channel.fromPath( params.comparison, checkIfExists: true ))
            )
            ch_dp = DIFFERENTIAL_PEAKS.out.data
            ch_versions = ch_versions.mix(DIFFERENTIAL_PEAKS.out.versions)
        }


    emit:
    versions = ch_versions
    dp = ch_dp

}
