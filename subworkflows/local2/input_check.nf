/*
 * Check input samplesheet and get read channels
 */

include { SAMPLESHEET_CHECK } from '../../modules/local2/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    metadata
    workflow

    main:
    SAMPLESHEET_CHECK (
        samplesheet,
        metadata,
        workflow
    )

    SAMPLESHEET_CHECK.out.csv
        .splitCsv ( header:true, sep:"," )
        .map { get_samplesheet_paths(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
    samplesheet = SAMPLESHEET_CHECK.out.csv
    versions = SAMPLESHEET_CHECK.out.versions
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
