/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local2/samplesheet_check'
include { GROUP_IGG } from '../../modules/local2/group_igg'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    metadata
    workflow
    use_igg

    main:
    ch_versions = Channel.empty()
    ch_fastq_multi = Channel.empty()
    ch_fastq_single = Channel.empty()
    meta_igg = Channel.empty()
    samplesheet_comb_igg = Channel.empty()

    SAMPLESHEET_CHECK (
        samplesheet,
        metadata,
        workflow
    )
    samplesheet = SAMPLESHEET_CHECK.out.csv

    fastq = samplesheet
        .splitCsv ( header:true, sep:"," )
        .map { get_samplesheet_paths(it) }
        .map {
            meta, fastq ->
                [ meta, fastq ]
            }
        .groupTuple(by: [0])

    ch_fastq_multi = fastq
        .filter { meta, fastq -> fastq.size() > 1 }
        .map { meta, fastq -> [ meta, fastq.flatten() ] }

    ch_fastq_single = fastq
        .filter { meta, fastq -> fastq.size() == 1 }
        .map { meta, fastq -> [ meta, fastq.flatten() ] }

    ch_versions = SAMPLESHEET_CHECK.out.versions


    /* process IgG */
    def use_igg_modes = use_igg instanceof Collection ? use_igg : [use_igg]
    def valid_igg_modes = ['group', 'all', 'custom']
    if (use_igg_modes.any { it in valid_igg_modes }) {
        def igg_group = use_igg_modes.find { it in valid_igg_modes }
        GROUP_IGG(
        samplesheet,
        igg_group
        )
        samplesheet_comb_igg = GROUP_IGG.out.csv

        GROUP_IGG.out.csv
            .splitCsv ( header:true, sep:"," )
            .map { format_csv(it) }
            .set { meta_igg }

        ch_versions = ch_versions.mix(GROUP_IGG.out.versions)
    }

    emit:
    fastq_multi = ch_fastq_multi // channel: [ val(meta), [ reads ] ]
    fastq_single = ch_fastq_single // channel: [ val(meta), [ reads ] ]
    samplesheet
    samplesheet_comb_igg
    meta_igg // channel: [ val(meta) ]
    versions = ch_versions

}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row) {
    def meta = row
    meta.single_end    = meta.single_end.toBoolean()
    meta.is_control    = meta.is_control.toBoolean()
    meta.call_peak       = meta.call_peak.toBoolean()
    meta.call_rep_peak   = meta.call_rep_peak.toBoolean()
    meta.call_con_peak   = meta.call_con_peak.toBoolean()


    def array = []
    if (!file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1) ] ]
    } else {
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        array = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    return array
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def format_csv(LinkedHashMap row) {
    def meta = row
    meta.single_end    = meta.single_end.toBoolean()
    meta.is_control    = meta.is_control.toBoolean()
    meta.call_peak       = meta.call_peak.toBoolean()
    meta.call_rep_peak   = meta.call_rep_peak.toBoolean()
    meta.call_con_peak   = meta.call_con_peak.toBoolean()

    return meta
}
