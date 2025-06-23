/*
 * Check input samplesheet and get read channels
 */

include { GROUP_IGG } from '../../modules/local2/group_igg'

workflow IGG_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    igg_group

    main:
    GROUP_IGG(
        samplesheet,
        igg_group
    )

    GROUP_IGG.out.csv
        .splitCsv ( header:true, sep:"," )
        .map { format_csv(it) }
        .set { meta }

    emit:
    meta // channel: [ val(meta) ]
    samplesheet = GROUP_IGG.out.csv
    versions = GROUP_IGG.out.versions
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
