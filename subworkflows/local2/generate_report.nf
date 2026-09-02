/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { MAKE_REPORT      } from '../../modules/local2/make_report'
include { MAKE_REPORT as MAKE_REPORT_COMB_IGG      } from '../../modules/local2/make_report'

workflow GENERATE_REPORT {
    take:
    ch_read_metrics
    ch_frag_lens
    samplesheet
    ch_orig_csv
    ch_orig_widths
    ch_rip
    ch_rep_csv
    ch_conp_csv
    ch_conp_bed
    ch_conp_ann
    ch_dp
    samplesheet_comb_igg
    ch_orig_csv_comb_igg
    ch_orig_widths_comb_igg
    ch_rip_comb_igg
    ch_rep_csv_comb_igg
    ch_conp_csv_comb_igg
    ch_conp_bed_comb_igg
    ch_conp_ann_comb_igg
    ch_dp_comb_igg
    use_igg
    report_dir
    srcdir


    main:
    ch_versions = Channel.empty()

    if (!params.run_read_qc){
        ch_read_metrics = Channel.fromPath("${srcdir}/cached_data/read_metrics.csv", checkIfExists: true)
        ch_frag_lens= Channel.fromPath("${srcdir}/cached_data/fragment_lengths/", type: "dir", checkIfExists: true)
    }

    if ('individual' in use_igg){
        if (!params.run_peak_qc){
            ch_orig_csv = Channel.fromPath("${srcdir}/cached_data/individual_igg/original_peaks/*", checkIfExists: true)
            ch_orig_widths = Channel.fromPath("${srcdir}/cached_data/individual_igg/replicated_peaks/*", checkIfExists: true)
            ch_rip = Channel.fromPath("${srcdir}/cached_data/individual_igg/reads_in_peak/*", checkIfExists: true)
            ch_rep_csv = Channel.fromPath("${srcdir}/cached_data/individual_igg/replicated_peaks/*", checkIfExists: true)
            ch_conp_csv = Channel.fromPath("${srcdir}/cached_data/individual_igg/consensus_peaks/*", checkIfExists: true)
            ch_conp_bed = Channel.fromPath("${srcdir}/csv/consensus_peaks.individual_igg.csv", checkIfExists: true)
                .splitCsv(header: true)
                .map { row -> [ row.target, file("${srcdir}/${row.conp_bed}", checkIfExists: true) ] }
            ch_conp_ann= Channel.fromPath("${srcdir}/2_consensus_peaks/01_peak_annotations/*", checkIfExists: true)
        }

        if (!params.run_differential_peaks){
            ch_dp = Channel.fromPath("${srcdir}/cached_data/individual_igg/differential_peaks/*", checkIfExists: true)
        }

        MAKE_REPORT(
                samplesheet,
                ch_read_metrics.collect().ifEmpty([]),
                ch_frag_lens.collect().ifEmpty([]),
                ch_orig_csv.collect().ifEmpty([]),
                ch_orig_widths.collect().ifEmpty([]),
                ch_rip.collect().ifEmpty([]),
                ch_rep_csv.collect().ifEmpty([]),
                ch_conp_csv.collect().ifEmpty([]),
                ch_conp_bed.collect{it[1]}.ifEmpty([]),
                ch_conp_ann.collect().ifEmpty([]),
                ch_dp.collect().ifEmpty([]),
                report_dir
        )
        ch_versions = MAKE_REPORT.out.versions
    }

    if ([ 'group', 'all', 'custom' ].any { it in use_igg }){
        if (!params.run_peak_qc){
            ch_orig_csv_comb_igg = Channel.fromPath("${srcdir}/cached_data/combined_igg/original_peaks/*", checkIfExists: true)
            ch_orig_widths_comb_igg = Channel.fromPath("${srcdir}/cached_data/combined_igg/replicated_peaks/*", checkIfExists: true)
            ch_rip_comb_igg = Channel.fromPath("${srcdir}/cached_data/combined_igg/reads_in_peak/*", checkIfExists: true)
            ch_rep_csv_comb_igg = Channel.fromPath("${srcdir}/cached_data/combined_igg/replicated_peaks/*", checkIfExists: true)
            ch_conp_csv_comb_igg = Channel.fromPath("${srcdir}/cached_data/combined_igg/consensus_peaks/*", checkIfExists: true)
            ch_conp_bed_comb_igg = Channel.fromPath("${srcdir}/csv/consensus_peaks.combined_igg.csv", checkIfExists: true)
                .splitCsv(header: true)
                .map { row -> [ row.target, file("${srcdir}/${row.conp_bed}", checkIfExists: true) ] }
            ch_conp_ann_comb_igg = Channel.fromPath("${srcdir}/2_consensus_peaks_combined_igg/01_peak_annotations/*", checkIfExists: true)
        }

        if (!params.run_differential_peaks){
            ch_dp_comb_igg = Channel.fromPath("${srcdir}/cached_data/combined_igg/differential_peaks/*", checkIfExists: true)
        }

        MAKE_REPORT_COMB_IGG(
                samplesheet_comb_igg,
                ch_read_metrics.collect().ifEmpty([]),
                ch_frag_lens.collect().ifEmpty([]),
                ch_orig_csv_comb_igg.collect().ifEmpty([]),
                ch_orig_widths_comb_igg.collect().ifEmpty([]),
                ch_rip_comb_igg.collect().ifEmpty([]),
                ch_rep_csv_comb_igg.collect().ifEmpty([]),
                ch_conp_csv_comb_igg.collect().ifEmpty([]),
                ch_conp_bed_comb_igg.collect{it[1]}.ifEmpty([]),
                ch_conp_ann_comb_igg.collect().ifEmpty([]),
                ch_dp_comb_igg.collect().ifEmpty([]),
                report_dir
        )
        ch_versions = MAKE_REPORT_COMB_IGG.out.versions
    }

    emit:
    versions = ch_versions
}
