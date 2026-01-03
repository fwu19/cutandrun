/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

// Validate input parameters in specialised library
WorkflowCutandrun.initialise(params, log)
def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

// Check input path parameters to see if the files exist if they have been specified
checkPathParamList = [
    params.fasta,
    params.gtf
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check spike-in
checkPathParamList = [
        params.spikein_bowtie2,
        params.spikein_fasta
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Check mandatory parameters that cannot be checked in the groovy lib as we want a channel for them
//if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified!" }

ch_blacklist = Channel.empty()
if (params.blacklist) {
    ch_blacklist = Channel.from( file(params.blacklist) )
}
else {
    ch_blacklist = Channel.empty()
    WorkflowCutandrun.blacklistWarn(log)
}

// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
ch_dummy_csv = file("$projectDir/assets/dummy_file.csv", checkIfExists: true)

// Stage awk files for parsing log files
ch_bt2_to_csv_awk     = file("$projectDir/bin/bt2_report_to_csv.awk"    , checkIfExists: true)
ch_dt_frag_to_csv_awk = file("$projectDir/bin/dt_frag_report_to_csv.awk", checkIfExists: true)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

// Load up and check multiqc base config and custom configs
ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.fromPath("$projectDir/assets/local/multiqc_config.yml")


/*
========================================================================================
    INIALISE PARAMETERS AND OPTIONS
========================================================================================
*/

// Init aligners
def prepare_tool_indices = ["bowtie2"]

// Check peak caller params
def caller_list = ['seacr', 'macs2']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : ['macs2', 'seacr']
if ((caller_list + callers).unique().size() != caller_list.size()) {
    exit 1, "Invalid variant calller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}"
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * SUBWORKFLOWS
 */
include { PREPARE_GENOME                                   } from "../subworkflows/local2/prepare_genome"
include { FASTQC_TRIMGALORE                                } from "../subworkflows/local/fastqc_trimgalore"
include { ALIGN_BOWTIE2                                    } from "../subworkflows/local/align_bowtie2"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_TARGET_META  } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_SPIKEIN_META } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_PICARD_DUP_META  } from "../subworkflows/local/extract_metadata_awk"
include { MARK_DUPLICATES_PICARD                           } from "../subworkflows/local/mark_duplicates_picard"
include { MARK_DUPLICATES_PICARD as DEDUPLICATE_PICARD     } from "../subworkflows/local/mark_duplicates_picard"
include { SAMTOOLS_VIEW_SORT_STATS as FILTER_READS         } from "../subworkflows/local/samtools_view_sort_stats"

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { CAT_FASTQ                                                    } from "../modules/nf-core/cat/fastq/main"

/*
========================================================================================
    IMPORT CUSTOM MODULES/SUBWORKFLOWS
========================================================================================
*/
include { INPUT_CHECK                                                        } from "../subworkflows/local2/input_check"
include { IGG_CHECK                                                          } from '../subworkflows/local2/igg_check'
include { CALL_SEACR_PEAKS                                                   } from '../subworkflows/local2/call_seacr_peaks'
include { CALL_MACS2_PEAKS                                                   } from '../subworkflows/local2/call_macs2_peaks'
include { CALL_PEAKS_PROCESS_CONTROLS                                        } from '../subworkflows/local2/call_peaks_process_controls'
include { COMPUTE_GENOMECOVERAGE                                             } from "../subworkflows/local2/compute_genomecoverage"
include { QC_READS                                                           } from '../subworkflows/local2/qc_reads'
include { QC_PEAKS                                                           } from '../subworkflows/local2/qc_peaks'
include { QC_PEAKS as QC_PEAKS_COMB_IGG                                      } from '../subworkflows/local2/qc_peaks'
include { QC_PROCESS_CONTROLS                                                } from '../subworkflows/local2/qc_process_controls'
include { SUMMARY_PLOTS                                                      } from '../subworkflows/local2/summary_plots'
include { SUMMARY_PLOTS as SUMMARY_PLOTS_COMB_IGG                            } from '../subworkflows/local2/summary_plots'

include { GET_FASTQ_PATHS                                                    } from '../modules/local2/get_fastq_paths'
include { MULTIQC                                                            } from '../modules/local2/multiqc'
include { READS_IN_CONSENSUS_PEAKS                                           } from '../modules/local2/reads_in_consensus_peaks'
include { READS_IN_CONSENSUS_PEAKS as READS_IN_CONSENSUS_PEAKS_COMB_IGG      } from '../modules/local2/reads_in_consensus_peaks'
include { COLLECT_COUNT_MATRIX                                               } from '../modules/local2/collect_count_matrix'
include { COLLECT_COUNT_MATRIX as COLLECT_COUNT_MATRIX_COMB_IGG              } from '../modules/local2/collect_count_matrix'
include { DIFFERENTIAL_PEAKS                                                 } from '../modules/local2/differential_peaks'
include { DIFFERENTIAL_PEAKS as DIFFERENTIAL_PEAKS_COMB_IGG                  } from '../modules/local2/differential_peaks'
include { GENERATE_REPORT                                                    } from '../modules/local2/generate_report'
include { GENERATE_REPORT as GENERATE_REPORT_COMB_IGG                        } from '../modules/local2/generate_report'
include { GENERATE_REPORT_PROCESS_CONTROLS                                   } from '../modules/local2/generate_report_process_controls'
include { WRITE_CSV as WRITE_CSV_BT2                                         } from '../modules/local2/write_csv'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CUTANDRUN {

    // Init
    ch_software_versions = Channel.empty()

    /*
     * SUBWORKFLOW: Uncompress and prepare reference genome files
     */
    if (!params.target_genome){ params.target_genome = params.genome }

    if(params.run_genome_prep) {
        PREPARE_GENOME (
            prepare_tool_indices,
            ch_blacklist
        )
        ch_software_versions = ch_software_versions.mix(PREPARE_GENOME.out.versions)
    }

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    if(params.run_input_check) {

        /* Get fastq paths */
        if ( params.input_dir =~ 'dummy' ){
            if ( params.input =~ 'dummy' ){
                exit 1, 'Neither --input nor --input_dir is specified!'
            }else {
                ch_input = Channel.fromPath( params.input, checkIfExists: true )
            }
        }else {
            GET_FASTQ_PATHS (
                Channel.fromPath("${params.input_dir}", checkIfExists: true),
                params.workflow
            )
            ch_input = GET_FASTQ_PATHS.out.csv
            ch_software_versions = ch_software_versions.mix(GET_FASTQ_PATHS.out.versions)
        }

        /* Add metadata */
        ch_metadata = params.metadata ? file( params.metadata, checkIfExists: true ) : ch_dummy_csv
        INPUT_CHECK (
            ch_input,
            ch_metadata,
            params.workflow
        )

        samplesheet = INPUT_CHECK.out.samplesheet
        ch_software_versions = ch_software_versions.mix(INPUT_CHECK.out.versions)

        /* Generate sample sheet */
        INPUT_CHECK.out.reads
        .map {
            meta, fastq ->
                [ meta, fastq ] }
        .groupTuple(by: [0])
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }

        /* Update sample sheet with combined IgG */
        samplesheet_combine_igg = Channel.empty()
        meta_combine_igg = Channel.empty()
        if (params.run_combine_igg){
            IGG_CHECK(
                samplesheet,
                params.igg_group
            )
            samplesheet_combine_igg = IGG_CHECK.out.samplesheet
            meta_combine_igg = IGG_CHECK.out.meta
            ch_software_versions = ch_software_versions.mix(IGG_CHECK.out.versions)
        }
    }


    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    if(params.run_cat_fastq) {
        CAT_FASTQ (
            ch_fastq.multiple
        )
        ch_software_versions = ch_software_versions.mix(CAT_FASTQ.out.versions)

        CAT_FASTQ.out.reads
        .mix(ch_fastq.single)
        .set { ch_cat_fastq }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [READS]]
    //ch_cat_fastq | view

    /*
     * SUBWORKFLOW: Read QC, trim adapters and perform post-trim read QC
     */
    if(params.run_trim_galore_fastqc) {
        FASTQC_TRIMGALORE (
            ch_cat_fastq,
            params.skip_fastqc,
            params.skip_trimming
        )
        ch_trimmed_reads     = FASTQC_TRIMGALORE.out.reads
        ch_software_versions = ch_software_versions.mix(FASTQC_TRIMGALORE.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [READS]]
    //FASTQC_TRIMGALORE.out.reads | view

    /*
    * SUBWORKFLOW: Alignment to target and spikein genome using botwtie2
    */
    ch_orig_bam                   = Channel.empty()
    ch_orig_spikein_bam           = Channel.empty()
    ch_bowtie2_log                = Channel.empty()
    ch_bowtie2_spikein_log        = Channel.empty()
    ch_samtools_bam               = Channel.empty()
    ch_samtools_bai               = Channel.empty()
    ch_samtools_stats             = Channel.empty()
    ch_samtools_flagstat          = Channel.empty()
    ch_samtools_idxstats          = Channel.empty()
    ch_samtools_spikein_bam       = Channel.empty()
    ch_samtools_spikein_bai       = Channel.empty()
    ch_samtools_spikein_stats     = Channel.empty()
    ch_samtools_spikein_flagstat  = Channel.empty()
    ch_samtools_spikein_idxstats  = Channel.empty()
    if(params.run_alignment) {
        if (params.aligner == "bowtie2") {
            ALIGN_BOWTIE2 (
                ch_trimmed_reads,
                PREPARE_GENOME.out.bowtie2_index,
                PREPARE_GENOME.out.bowtie2_spikein_index,
                PREPARE_GENOME.out.fasta,
                PREPARE_GENOME.out.spikein_fasta
            )
            ch_orig_bam                   = ALIGN_BOWTIE2.out.orig_bam
            ch_orig_spikein_bam           = ALIGN_BOWTIE2.out.orig_spikein_bam
            ch_bowtie2_log                = ALIGN_BOWTIE2.out.bowtie2_log
            ch_bowtie2_spikein_log        = ALIGN_BOWTIE2.out.bowtie2_spikein_log

            ch_samtools_bam               = ALIGN_BOWTIE2.out.bam
            ch_samtools_bai               = ALIGN_BOWTIE2.out.bai
            ch_samtools_stats             = ALIGN_BOWTIE2.out.stats
            ch_samtools_flagstat          = ALIGN_BOWTIE2.out.flagstat
            ch_samtools_idxstats          = ALIGN_BOWTIE2.out.idxstats

            ch_samtools_spikein_bam       = ALIGN_BOWTIE2.out.spikein_bam
            ch_samtools_spikein_bai       = ALIGN_BOWTIE2.out.spikein_bai
            ch_samtools_spikein_stats     = ALIGN_BOWTIE2.out.spikein_stats
            ch_samtools_spikein_flagstat  = ALIGN_BOWTIE2.out.spikein_flagstat
            ch_samtools_spikein_idxstats  = ALIGN_BOWTIE2.out.spikein_idxstats
            ch_software_versions          = ch_software_versions.mix(ALIGN_BOWTIE2.out.versions)

        }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: extract aligner metadata
     */
    ch_metadata_bt2_target  = Channel.empty()
    ch_metadata_bt2_spikein = Channel.empty()
    if (params.aligner == "bowtie2" && params.run_alignment) {
        EXTRACT_BT2_TARGET_META (
            ch_bowtie2_log,
            ch_bt2_to_csv_awk,
            true
        )
        ch_metadata_bt2_target = EXTRACT_BT2_TARGET_META.out.metadata
        ch_software_versions   = ch_software_versions.mix(EXTRACT_BT2_TARGET_META.out.versions)

        EXTRACT_BT2_SPIKEIN_META (
            ch_bowtie2_spikein_log,
            ch_bt2_to_csv_awk,
            true
        )
        ch_metadata_bt2_spikein = EXTRACT_BT2_SPIKEIN_META.out.metadata
    }
    //ch_metadata_bt2_target | view
    //ch_metadata_bt2_spikein | view

    /*
     *  SUBWORKFLOW: Filter reads based some standard measures
     *  - Unmapped reads 0x004
     *  - Mate unmapped 0x0008
     *  - Multi-mapped reads
     *  - Filter out reads aligned to blacklist regions
     *  - Filter out reads below a threshold q score
     *  - Filter out mitochondrial reads (if required)
     */
    if (params.run_read_filter) {
        FILTER_READS (
            ch_samtools_bam,
            PREPARE_GENOME.out.allowed_regions.collect{it[1]}.ifEmpty([]),
            PREPARE_GENOME.out.fasta
        )
        ch_samtools_bam      = FILTER_READS.out.bam
        ch_samtools_bai      = FILTER_READS.out.bai
        ch_samtools_stats    = FILTER_READS.out.stats
        ch_samtools_flagstat = FILTER_READS.out.flagstat
        ch_samtools_idxstats = FILTER_READS.out.idxstats
        ch_software_versions = ch_software_versions.mix(FILTER_READS.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: Mark duplicates on all samples
     */
    ch_samtools_bam_markdup = Channel.empty()
    ch_samtools_bai_markdup = Channel.empty()
    ch_samtools_bam_dedup = Channel.empty()
    ch_samtools_bai_dedup = Channel.empty()
    ch_markduplicates_metrics = Channel.empty()

    if (params.run_mark_dups) {
        MARK_DUPLICATES_PICARD (
            ch_samtools_bam,
            ch_samtools_bai,
            true,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect()
        )
        ch_samtools_bam_markdup           = MARK_DUPLICATES_PICARD.out.bam
        ch_samtools_bai_markdup           = MARK_DUPLICATES_PICARD.out.bai
        ch_samtools_stats_markdup         = MARK_DUPLICATES_PICARD.out.stats
        ch_samtools_flagstat_markdup      = MARK_DUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats_markdup      = MARK_DUPLICATES_PICARD.out.idxstats
        ch_markduplicates_metrics = MARK_DUPLICATES_PICARD.out.metrics
        ch_software_versions      = ch_software_versions.mix(MARK_DUPLICATES_PICARD.out.versions)

        // write out metadata with paths to bam files
        WRITE_CSV_BT2(
            ch_samtools_bam_markdup
                .join(ch_samtools_bai_markdup)
                .map {
                    meta, bam, bai -> meta + [target_bam: "${params.outdir}/1_individual_samples/02_alignment/bowtie2/${params.target_genome}/${bam.name}"] + [target_bai: "${params.outdir}/1_individual_samples/02_alignment/bowtie2/${params.target_genome}/${bai.name}"]
                    }
                .collect(),
            "bt2_bam_markdup.csv"
        )

    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bam | view

    /*
     * SUBWORKFLOW: Remove duplicates - default is on IgG controls only
     */
    if (params.run_remove_dups) {
        DEDUPLICATE_PICARD (
            ch_samtools_bam,
            ch_samtools_bai,
            params.dedup_target_reads,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect()
        )
        ch_samtools_bam_dedup      = DEDUPLICATE_PICARD.out.bam
        ch_samtools_bai_dedup      = DEDUPLICATE_PICARD.out.bai
        ch_samtools_stats_dedup    = DEDUPLICATE_PICARD.out.stats
        ch_samtools_flagstat_dedup = DEDUPLICATE_PICARD.out.flagstat
        ch_samtools_idxstats_dedup = DEDUPLICATE_PICARD.out.idxstats
        ch_software_versions = ch_software_versions.mix(DEDUPLICATE_PICARD.out.versions)
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [BAM]]
    //ch_samtools_bai | view

    /*
    * SUBWORKFLOW: extract duplication stats from picard report
    */
    ch_metadata_picard_duplicates = Channel.empty()
    if (params.run_mark_dups) {
        EXTRACT_PICARD_DUP_META (
            ch_markduplicates_metrics,
            ch_dummy_file.collect(),
            false
        )
        ch_metadata_picard_duplicates = EXTRACT_PICARD_DUP_META.out.metadata
        ch_software_versions          = ch_software_versions.mix(EXTRACT_PICARD_DUP_META.out.versions)
    }
    //ch_metadata_picard_duplicates | view


    /*
    * SUBWORKFLOW: Convert BAM files to bedgraph/bigwig and apply spikein normalisation if required
    */
    ch_bedgraph_markdup     = Channel.empty()
    ch_bedgraph_dedup       = Channel.empty()
    ch_bigwig_markdup       = Channel.empty()
    ch_bigwig_dedup         = Channel.empty()

    if(params.run_mark_dups && params.run_remove_dups) {
        COMPUTE_GENOMECOVERAGE(
            ch_samtools_bam_markdup,
            ch_samtools_bai_markdup,
            ch_samtools_bam_dedup,
            ch_samtools_bai_dedup,
            PREPARE_GENOME.out.chrom_sizes.collect(),
            ch_dummy_file,
            params.normalisation_mode,
            ch_metadata_bt2_spikein
        )
        ch_bedgraph_markdup          = COMPUTE_GENOMECOVERAGE.out.bedgraph_markdup_unnorm
        ch_bedgraph_dedup          = COMPUTE_GENOMECOVERAGE.out.bedgraph_dedup_unnorm
        ch_bigwig_markdup            = COMPUTE_GENOMECOVERAGE.out.bigwig_markdup
        ch_bigwig_dedup            = COMPUTE_GENOMECOVERAGE.out.bigwig_dedup
        ch_software_versions = ch_software_versions.mix(COMPUTE_GENOMECOVERAGE.out.versions)

    }

    /*
    * Run MultiQC
    */
    ch_multiqc_data = Channel.empty()
    if (params.run_multiqc){

        MULTIQC(
            ch_multiqc_custom_config.ifEmpty([]),
            ch_bowtie2_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_spikein_log.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_metrics.collect{it[1]}.ifEmpty([])

        )
        ch_multiqc_data = MULTIQC.out.data
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions)
    }

    /*
     * SUBWORKFLOW: Call peaks from individual samples
     */
    if(params.run_peak_calling & params.workflow == "cutandrun") {

        CALL_SEACR_PEAKS (
            ch_bedgraph_markdup,
            ch_bedgraph_dedup,
            params.run_combine_igg,
            meta_combine_igg,
            params.skip_individual_igg
        )

        ch_software_versions = ch_software_versions.mix(CALL_SEACR_PEAKS.out.versions)

        CALL_MACS2_PEAKS (
            ch_samtools_bam_markdup,
            params.run_combine_igg,
            meta_combine_igg,
            params.skip_individual_igg
        )
        ch_software_versions = ch_software_versions.mix(CALL_MACS2_PEAKS.out.versions)

    }

    if(params.run_peak_calling & params.workflow == "process_controls") {
        CALL_PEAKS_PROCESS_CONTROLS (
            ch_bedgraph_markdup,
            ch_bedgraph_dedup,
            ch_samtools_bam_markdup,
            params.igg_dir,
            params.use_igg
        )
        ch_peaks_all = CALL_PEAKS_PROCESS_CONTROLS.out.peaks_all
        // ch_peaks_all.view()
        // [ meta, [peaks] ]

        ch_peaks_final = CALL_PEAKS_PROCESS_CONTROLS.out.peaks_final
        // ch_peaks_final.view()
        // [ meta, [peaks] ]

        ch_software_versions = ch_software_versions.mix(CALL_PEAKS_PROCESS_CONTROLS.out.versions)
    }

    if (params.run_read_qc){
        QC_READS(
            ch_multiqc_data,
            ch_samtools_bam
        )
        ch_software_versions = ch_software_versions.mix(QC_READS.out.versions)
    }

    /*
    * QC peaks
    */
    if (params.run_peak_qc && params.workflow == "cutandrun" && !params.skip_individual_igg){
        CALL_MACS2_PEAKS.out.narrow_filtered
            .concat(
                CALL_MACS2_PEAKS.out.narrow_igg,
                CALL_MACS2_PEAKS.out.narrow_noigg,
                CALL_MACS2_PEAKS.out.broad_filtered,
                CALL_MACS2_PEAKS.out.broad_igg,
                CALL_MACS2_PEAKS.out.broad_noigg,
                CALL_SEACR_PEAKS.out.seacr_filtered,
                CALL_SEACR_PEAKS.out.seacr_igg,
                CALL_SEACR_PEAKS.out.seacr_noigg
        )
        .groupTuple(by: 0)
        .set { ch_peaks_all}
        // ch_peaks_all.view()
        // [ [meta], [path(peak1), path(peak2), ...] ]


        CALL_MACS2_PEAKS.out.narrow_filtered
            .mix( CALL_MACS2_PEAKS.out.narrow_noigg.filter(it -> it[0].control_group == "") )
            .concat (
            CALL_MACS2_PEAKS.out.broad_filtered
                .mix( CALL_MACS2_PEAKS.out.broad_noigg.filter(it -> it[0].control_group == "") ),
            CALL_SEACR_PEAKS.out.seacr_filtered
                .mix( CALL_SEACR_PEAKS.out.seacr_noigg.filter(it -> it[0].control_group == "") )
            )
            .groupTuple(by: 0)
            .set { ch_peaks_final}
        // ch_peaks_final.view()
        // [ [meta], path(macs2_narrow_peak), path(macs2_broad_peak), path(seacr_peak) ]

        QC_PEAKS(
            samplesheet,
            params.min_replicates,
            params.fasta,
            params.gtf ? file(params.gtf, checkIfExists: true) : ch_dummy_file,
            ch_samtools_bam,
            ch_peaks_all,
            ch_peaks_final

        )
        ch_software_versions = ch_software_versions.mix(QC_PEAKS.out.versions)
    }

    /*
    * QC peaks with combined IgG
    */
    if (params.run_peak_qc && params.workflow == "cutandrun" && params.run_combine_igg){
        CALL_MACS2_PEAKS.out.narrow_comb_igg_filtered
            .concat(
                CALL_MACS2_PEAKS.out.narrow_comb_igg,
                CALL_MACS2_PEAKS.out.narrow_noigg,
                CALL_MACS2_PEAKS.out.broad_comb_igg_filtered,
                CALL_MACS2_PEAKS.out.broad_comb_igg,
                CALL_MACS2_PEAKS.out.broad_noigg,
                CALL_SEACR_PEAKS.out.seacr_comb_igg_filtered,
                CALL_SEACR_PEAKS.out.seacr_comb_igg,
                CALL_SEACR_PEAKS.out.seacr_noigg
        )
        .groupTuple(by: 0)
        .set { ch_peaks_comb_igg_all}
        // ch_peaks_all.view()
        // [ [meta], [path(peak1), path(peak2), ...] ]


        CALL_MACS2_PEAKS.out.narrow_comb_igg_filtered
            .mix( CALL_MACS2_PEAKS.out.narrow_noigg.filter(it -> it[0].control_group == "") )
            .concat (
            CALL_MACS2_PEAKS.out.broad_comb_igg_filtered
                .mix( CALL_MACS2_PEAKS.out.broad_noigg.filter(it -> it[0].control_group == "") ),
            CALL_SEACR_PEAKS.out.seacr_comb_igg_filtered
                .mix( CALL_SEACR_PEAKS.out.seacr_noigg.filter(it -> it[0].control_group == "") )
            )
            .groupTuple(by: 0)
            .set { ch_peaks_comb_igg_final}
        // ch_peaks_final.view()
        // [ [meta], path(macs2_narrow_peak), path(macs2_broad_peak), path(seacr_peak) ]

        QC_PEAKS_COMB_IGG(
            samplesheet_combine_igg,
            params.min_replicates,
            params.fasta,
            params.gtf ? file(params.gtf, checkIfExists: true) : ch_dummy_file,
            ch_samtools_bam,
            ch_peaks_comb_igg_all,
            ch_peaks_comb_igg_final

        )
        ch_software_versions = ch_software_versions.mix(QC_PEAKS_COMB_IGG.out.versions)
    }

    /*
    * QC peaks for process controls
    */
    if (params.run_peak_qc && params.workflow == "process_controls"){
        ch_ref_peaks = params.ref_peaks ? Channel.fromPath("${params.ref_peaks}", type: "dir", checkIfExists: true) : Channel.empty()
        QC_PROCESS_CONTROLS(
            samplesheet,
            params.genome,
            ch_samtools_bam,
            ch_peaks_all,
            ch_peaks_final,
            ch_ref_peaks.ifEmpty([])
        )
        ch_software_versions = ch_software_versions.mix(QC_PROCESS_CONTROLS.out.versions)
    }

    /*
    * Quantify reads in consensus peaks and call differential peaks
    */
    ch_dp = Channel.empty()
    if (params.run_peak_qc && !params.skip_individual_igg){
        // Count reads in consensus peaks
        READS_IN_CONSENSUS_PEAKS(
            QC_PEAKS.out.conp_bed
                .cross (
                        ch_samtools_bam
                            .filter { it[0].target != "IgG" }
                            .map { it -> [ it[0].target, it ] }
                )
                .map { it -> [ it[1][1][0], it[1][1][1], it[0][1] ] }
        )
        ch_software_versions = ch_software_versions.mix(READS_IN_CONSENSUS_PEAKS.out.versions)

        ch_cts = Channel.empty()
        COLLECT_COUNT_MATRIX(
            READS_IN_CONSENSUS_PEAKS.out.count
                .groupTuple( by: 0 )
                .map { it -> [ it[0], it[1].flatten().collect() ]}
                .cross ( QC_PEAKS.out.conp_bed )
                .map { it -> [ it[0][0], it[0][1], it[1][1] ]}
                .combine(samplesheet)
        )
        ch_cts = COLLECT_COUNT_MATRIX.out.cts
        ch_software_versions = ch_software_versions.mix(COLLECT_COUNT_MATRIX.out.versions)

        // if --comparison is a dummy file or empty file is used, write an error message and continue.
        // if file has contents but not a correct format, throw an error.
        if ( params.run_differential_peaks ){
            DIFFERENTIAL_PEAKS(
                ch_cts
                .combine(samplesheet)
                .combine(Channel.fromPath( params.comparison, checkIfExists: true ))
            )
            ch_dp = DIFFERENTIAL_PEAKS.out.data
            ch_software_versions = ch_software_versions.mix(DIFFERENTIAL_PEAKS.out.versions)
        }
    }

    /*
    * Call differential peaks with combined IgG
    */
    ch_dp_comb_igg = Channel.empty()
    if (params.run_peak_qc && params.run_combine_igg){
        // Count reads in consensus peaks
        READS_IN_CONSENSUS_PEAKS_COMB_IGG(
            QC_PEAKS_COMB_IGG.out.conp_bed
                .cross (
                        ch_samtools_bam
                            .filter { it[0].target != "IgG" }
                            .map { it -> [ it[0].target, it ] }
                )
                .map { it -> [ it[1][1][0], it[1][1][1], it[0][1] ] }
        )
        ch_software_versions = ch_software_versions.mix(QC_PEAKS_COMB_IGG.out.versions)

        ch_cts_comb_igg = Channel.empty()
        COLLECT_COUNT_MATRIX_COMB_IGG(
            READS_IN_CONSENSUS_PEAKS_COMB_IGG.out.count
                .groupTuple( by: 0 )
                .map { it -> [ it[0], it[1].flatten().collect() ]}
                .cross ( QC_PEAKS_COMB_IGG.out.conp_bed )
                .map { it -> [ it[0][0], it[0][1], it[1][1] ]}
                .combine(samplesheet)
        )
        ch_cts_comb_igg = COLLECT_COUNT_MATRIX_COMB_IGG.out.cts
        ch_software_versions = ch_software_versions.mix(COLLECT_COUNT_MATRIX_COMB_IGG.out.versions)

        // if --comparison is a dummy file or empty file is used, throw a warning and continue.
        // if file has contents but not a correct format, throw an error.
        if( params.run_differential_peaks){
            DIFFERENTIAL_PEAKS_COMB_IGG(
                ch_cts_comb_igg
                .combine(samplesheet)
                .combine(Channel.fromPath( params.comparison, checkIfExists: true ))
            )
            ch_dp_comb_igg = DIFFERENTIAL_PEAKS_COMB_IGG.out.data
            ch_software_versions = ch_software_versions.mix(DIFFERENTIAL_PEAKS_COMB_IGG.out.versions)
        }
    }

    /*
    * make plots, e.g. heatmaps
    */
    if (params.run_peak_qc && params.run_summary_plots && !params.skip_individual_igg){
        SUMMARY_PLOTS(
            ch_bigwig_markdup,
            QC_PEAKS.out.conp_bed.ifEmpty([]),
            QC_PEAKS.out.conp_ann.collect{it[1]}.ifEmpty([]),
            ch_dp.collect().ifEmpty([]),
            params.gtf ? file(params.gtf, checkIfExists: true) : ch_dummy_file,
            params.gene_bed ? file(params.gene_bed, checkIfExists: true) : "$projectDir/assets/dummy_file.txt",
            true
        )
        ch_software_versions = ch_software_versions.mix(SUMMARY_PLOTS.out.versions)
    }

    if (params.run_peak_qc && params.run_summary_plots && params.run_combine_igg){
        SUMMARY_PLOTS_COMB_IGG(
            ch_bigwig_markdup,
            QC_PEAKS_COMB_IGG.out.conp_bed.ifEmpty([]),
            QC_PEAKS_COMB_IGG.out.conp_ann.collect{it[1]}.ifEmpty([]),
            ch_dp_comb_igg.collect().ifEmpty([]),
            params.gtf ? file(params.gtf, checkIfExists: true) : ch_dummy_file,
            params.gene_bed ? file(params.gene_bed, checkIfExists: true) : ch_dummy_file,
            false
        )
        ch_software_versions = ch_software_versions.mix(SUMMARY_PLOTS_COMB_IGG.out.versions)
    }

    /*
    * Generate report for experimental data
    */
    if (params.run_reporting && params.workflow == "cutandrun" && !params.skip_individual_igg){
            /*
            * Make plots for report
            */

            GENERATE_REPORT(
                samplesheet,
                QC_READS.out.read_metrics.ifEmpty([]),
                QC_READS.out.frag_lens.collect{it[1]}.ifEmpty([]),
                QC_PEAKS.out.orig_csv.collect{it[1]}.ifEmpty([]),
                QC_PEAKS.out.orig_widths.collect{it[1]}.ifEmpty([]),
                QC_PEAKS.out.rip.collect{it[1]}.ifEmpty([]),
                QC_PEAKS.out.rep_csv.collect{it[1]}.ifEmpty([]),
                QC_PEAKS.out.conp_csv.collect{it[1]}.ifEmpty([]),
                QC_PEAKS.out.conp_bed.collect{it[1]}.flatten().collect().ifEmpty([]),
                QC_PEAKS.out.conp_ann.collect{it[1]}.ifEmpty([]),
                ch_dp.flatten().collect().ifEmpty([]),
                params.report_dir ? Channel.fromPath("${params.report_dir}", type: 'dir', checkIfExists: true) : Channel.fromPath("$projectDir/assets/local/report/", type: 'dir', checkIfExists: true)
            )
            ch_software_versions = ch_software_versions.mix(GENERATE_REPORT.out.versions)

    }

    /*
    * Generate report for experimental data with combined IgG
    */
    if (params.run_reporting && params.workflow == "cutandrun" && params.run_combine_igg){
            /*
            * Make plots for report
            */

            GENERATE_REPORT_COMB_IGG(
                samplesheet_combine_igg,
                QC_READS.out.read_metrics.ifEmpty([]),
                QC_READS.out.frag_lens.collect{it[1]}.ifEmpty([]),
                QC_PEAKS_COMB_IGG.out.orig_csv.collect{it[1]}.ifEmpty([]),
                QC_PEAKS_COMB_IGG.out.orig_widths.collect{it[1]}.ifEmpty([]),
                QC_PEAKS_COMB_IGG.out.rip.collect{it[1]}.ifEmpty([]),
                QC_PEAKS_COMB_IGG.out.rep_csv.collect{it[1]}.ifEmpty([]),
                QC_PEAKS_COMB_IGG.out.conp_csv.collect{it[1]}.ifEmpty([]),
                QC_PEAKS_COMB_IGG.out.conp_bed.collect{it[1]}.flatten().collect().ifEmpty([]),
                QC_PEAKS_COMB_IGG.out.conp_ann.collect{it[1]}.ifEmpty([]),
                ch_dp_comb_igg.flatten().collect().ifEmpty([]),
                params.report_dir ? Channel.fromPath("${params.report_dir}", type: 'dir', checkIfExists: true) : Channel.fromPath("$projectDir/assets/local/report/", type: 'dir', checkIfExists: true)
            )
            ch_software_versions = ch_software_versions.mix(GENERATE_REPORT_COMB_IGG.out.versions)


    }

    /*
    * Generate report for process controls
    */
    if (params.run_reporting && params.workflow == "process_controls"){

            GENERATE_REPORT_PROCESS_CONTROLS(
                samplesheet,
                QC_READS.out.read_metrics.ifEmpty([]),
                QC_READS.out.frag_lens.collect{it[1]}.ifEmpty([]),
                QC_PROCESS_CONTROLS.out.orig_csv.collect{it[1]}.ifEmpty([]),
                QC_PROCESS_CONTROLS.out.orig_widths.collect{it[1]}.ifEmpty([]),
                QC_PROCESS_CONTROLS.out.rip.collect{it[1]}.ifEmpty([]),
                QC_PROCESS_CONTROLS.out.rep_csv.collect{it[1]}.ifEmpty([]),
                params.saved_data ? Channel.fromPath("${params.saved_data}", type: "dir", checkIfExists: true) : Channel.empty(),
                params.report_rmd ? Channel.fromPath("${params.report_rmd}", checkIfExists: true) : Channel.empty()
            )
            ch_software_versions = ch_software_versions.mix(GENERATE_REPORT_PROCESS_CONTROLS.out.versions)

    }

    /*
    * Collect software versions
    */
    ch_software_versions
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'software_versions.yml', sort: true, newLine: true)


}



////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////
import groovy.json.JsonOutput
workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
    def jsonStr = JsonOutput.toJson(params)
    def pretty  = JsonOutput.prettyPrint(jsonStr)
    file("${params.outdir ?: '.'}/pipeline_info/params.json").text = pretty

}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}


////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
