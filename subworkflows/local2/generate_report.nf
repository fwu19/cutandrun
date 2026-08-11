/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { MAKE_REPORT      } from '../../modules/local2/make_report'

workflow GENERATE_REPORT {
    take:
    samplesheet
    ch_read_metrics
    ch_frag_lens
    ch_orig_csv
    ch_orig_widths
    ch_rip
    ch_rep_csv
    ch_conp_csv
    ch_conp_bed
    ch_conp_ann
    ch_dp
    report_dir
    srcdir


    main:
    /*
        ch_orig_csv = Channel.empty()
        ch_orig_widths = Channel.empty()
        ch_rip = Channel.empty()
        ch_rep_bed = Channel.empty()
        ch_rep_csv = Channel.empty()
        ch_conp_bed = Channel.empty()
        ch_conp_csv = Channel.empty()
        ch_conp_ann = Channel.empty()
    */
        ch_versions = Channel.empty()

        MAKE_REPORT(
                samplesheet,
                ch_read_metrics.ifEmpty([]),
                ch_frag_lens.collect{it[1]}.ifEmpty([]),
                ch_orig_csv.collect{it[1]}.ifEmpty([]),
                ch_orig_widths.collect{it[1]}.ifEmpty([]),
                ch_rip.collect{it[1]}.ifEmpty([]),
                ch_rep_csv.collect{it[1]}.ifEmpty([]),
                ch_conp_csv.collect{it[1]}.ifEmpty([]),
                ch_conp_bed.collect{it[1]}.flatten().collect().ifEmpty([]),
                ch_conp_ann.collect{it[1]}.ifEmpty([]),
                ch_dp.flatten().collect().ifEmpty([]),
                report_dir
        )
        ch_versions = MAKE_REPORT.out.versions



    emit:
    versions = ch_versions
}
