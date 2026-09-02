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

// Check IgG usage
def igg_list = ['individual', 'all', 'custom', 'saved', 'none']
use_igg = params.use_igg ? params.use_igg.split(',').collect{ it.trim().toLowerCase()} : ['individual']
if ((igg_list + use_igg).unique().size() != igg_list.size()) {
    exit 1, "Invalid IgG usage option: ${params.use_igg}. Valid options: ${igg_list.join(', ')}"
}

if (params.step && !['mapping', 'peak_calling', 'peak_qc', 'consensus_peaks', 'differential_peaks'].contains(params.step)) {
    log.error "Invalid step specified: ${params.step}. Must be one of: mapping, peak_calling, peak_qc, consensus_peaks, differential_peaks"
    System.exit(1)
}

if (params.update && !['mapping', 'peak_calling', 'peak_qc', 'consensus_peaks', 'differential_peaks', 'summary_plots', 'report'].contains(params.update)) {
    log.error "Invalid update step specified: ${params.update}. Must be one of: mapping, peak_calling, peak_qc, consensus_peaks, differential_peaks, summary_plots, report"
    System.exit(1)
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
include { CALL_SEACR_PEAKS                                                   } from '../subworkflows/local2/call_seacr_peaks'
include { CALL_MACS2_PEAKS                                                   } from '../subworkflows/local2/call_macs2_peaks'
include { CALL_PEAKS_PROCESS_CONTROLS                                        } from '../subworkflows/local2/call_peaks_process_controls'
include { COMPUTE_GENOMECOVERAGE                                             } from "../subworkflows/local2/compute_genomecoverage"
include { QC_READS                                                           } from '../subworkflows/local2/qc_reads'
include { QC_PEAKS                                                           } from '../subworkflows/local2/qc_peaks'
include { QC_PROCESS_CONTROLS                                                } from '../subworkflows/local2/qc_process_controls'
include { CALL_DIFFERENTIAL_PEAKS                                            } from '../subworkflows/local2/call_differential_peaks'
include { SUMMARY_PLOTS                                                      } from '../subworkflows/local2/summary_plots'
include { GENERATE_REPORT                                                    } from '../subworkflows/local2/generate_report'
include { WRITE_OUTPUT_CSV                                                   } from "../subworkflows/local2/write_output_csv"

include { GET_FASTQ_PATHS                                                    } from '../modules/local2/get_fastq_paths'
include { MULTIQC                                                            } from '../modules/local2/multiqc'
include { GENERATE_REPORT_PROCESS_CONTROLS                                   } from '../modules/local2/generate_report_process_controls'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow CUTANDRUN {

    def outdir  = new File("${params.outdir}").absolutePath
    def srcdir  = params.srcdir ? new File("${params.srcdir}").absolutePath : outdir

    /*
    * Write software versions to a yaml file and overridden params to a json file
    */
    //WRITE_PARAMS()

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
    samplesheet = Channel.empty()
    ch_fastq_multi = Channel.empty()
    ch_fastq_single = Channel.empty()
    samplesheet = Channel.empty()
    samplesheet_combine_igg = Channel.empty()
    meta_combine_igg = Channel.empty()
    if(params.run_input_check) {

        /* Get fastq paths */
        if ( params.input_dir =~ 'dummy' ){
            if ( params.input =~ 'dummy' ){
                exit 1, 'Neither --input nor --input_dir is specified!'
            }else {
                samplesheet = Channel.fromPath( params.input, checkIfExists: true )
            }
        }else {
            GET_FASTQ_PATHS (
                Channel.fromPath("${params.input_dir}", checkIfExists: true),
                params.workflow
            )
            samplesheet = GET_FASTQ_PATHS.out.csv
            ch_software_versions = ch_software_versions.mix(GET_FASTQ_PATHS.out.versions)
        }

        /* Add metadata and create fastq channels */
        ch_metadata = params.metadata ? file( params.metadata, checkIfExists: true ) : ch_dummy_csv
        INPUT_CHECK (
            samplesheet,
            ch_metadata,
            params.workflow,
            use_igg
        )
        ch_fastq_multi = INPUT_CHECK.out.fastq_multi
        ch_fastq_single = INPUT_CHECK.out.fastq_single
        samplesheet = INPUT_CHECK.out.samplesheet
        samplesheet_combine_igg = INPUT_CHECK.out.samplesheet_comb_igg
        meta_combine_igg = INPUT_CHECK.out.meta_igg
        ch_software_versions = ch_software_versions.mix(INPUT_CHECK.out.versions)

    }


    /*
     * MODULE: Concatenate FastQ files from same sample if required
     */
    if(params.run_cat_fastq) {
        CAT_FASTQ (
            ch_fastq_multi
        )
        ch_software_versions = ch_software_versions.mix(CAT_FASTQ.out.versions)

        CAT_FASTQ.out.reads
        .mix(ch_fastq_single)
        .set { ch_cat_fastq }
    }
    //EXAMPLE CHANNEL STRUCT: [[id:h3k27me3_R1, group:h3k27me3, replicate:1, single_end:false, is_control:false], [READS]]
    //ch_cat_fastq | view

    /*
     * SUBWORKFLOW: Read QC, trim adapters and perform post-trim read QC
     */
    ch_trimmed_reads = Channel.empty()
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
    * SUBWORKFLOW: Run MultiQC
    */
    ch_multiqc_data = Channel.empty()
    if (params.run_multiqc){

        MULTIQC(
            ch_multiqc_custom_config.ifEmpty([]),
            ch_bowtie2_log.collect{it[1]},
            ch_bowtie2_spikein_log.collect{it[1]},
            ch_samtools_stats.collect{it[1]},
            ch_samtools_flagstat.collect{it[1]},
            ch_samtools_idxstats.collect{it[1]},
            ch_markduplicates_metrics.collect{it[1]}

        )
        ch_multiqc_data = MULTIQC.out.data
        ch_software_versions = ch_software_versions.mix(MULTIQC.out.versions)
    }

    /*
    * SUBWORKFLOW: QC reads
    */
    ch_read_metrics = Channel.empty()
    ch_frag_lens = Channel.empty()
    if (params.run_read_qc){
        QC_READS(
            ch_multiqc_data,
            ch_samtools_bam
        )
        ch_read_metrics = QC_READS.out.read_metrics
        ch_frag_lens = QC_READS.out.frag_lens
        ch_software_versions = ch_software_versions.mix(QC_READS.out.versions)
    }

    /*
     * SUBWORKFLOW: Call peaks from individual samples
     */
    peaks_narrow_igg = Channel.empty()
    peaks_narrow_noigg = Channel.empty()
    peaks_narrow_filtered = Channel.empty()
    peaks_narrow_comb_igg_filtered = Channel.empty()
    peaks_narrow_comb_igg = Channel.empty()
    peaks_broad_igg = Channel.empty()
    peaks_broad_noigg = Channel.empty()
    peaks_broad_filtered = Channel.empty()
    peaks_broad_comb_igg_filtered = Channel.empty()
    peaks_broad_comb_igg = Channel.empty()
    peaks_seacr_igg = Channel.empty()
    peaks_seacr_noigg = Channel.empty()
    peaks_seacr_filtered = Channel.empty()
    peaks_seacr_comb_igg_filtered = Channel.empty()
    peaks_seacr_comb_igg = Channel.empty()
    if(params.run_peak_calling && params.workflow == "cutandrun") {
        CALL_SEACR_PEAKS (
            ch_bedgraph_markdup,
            ch_bedgraph_dedup,
            meta_combine_igg,
            use_igg,
            srcdir
        )
        peaks_seacr_igg = CALL_SEACR_PEAKS.out.seacr_igg
        peaks_seacr_noigg = CALL_SEACR_PEAKS.out.seacr_noigg
        peaks_seacr_filtered = CALL_SEACR_PEAKS.out.seacr_filtered
        peaks_seacr_comb_igg_filtered = CALL_SEACR_PEAKS.out.seacr_comb_igg_filtered
        peaks_seacr_comb_igg = CALL_SEACR_PEAKS.out.seacr_comb_igg
        ch_software_versions = ch_software_versions.mix(CALL_SEACR_PEAKS.out.versions)

        CALL_MACS2_PEAKS (
            ch_samtools_bam_markdup,
            meta_combine_igg,
            use_igg,
            srcdir
        )
        peaks_narrow_igg = CALL_MACS2_PEAKS.out.narrow_igg
        peaks_narrow_noigg = CALL_MACS2_PEAKS.out.narrow_noigg
        peaks_narrow_filtered = CALL_MACS2_PEAKS.out.narrow_filtered
        peaks_narrow_comb_igg_filtered = CALL_MACS2_PEAKS.out.narrow_comb_igg_filtered
        peaks_narrow_comb_igg = CALL_MACS2_PEAKS.out.narrow_comb_igg
        peaks_broad_igg = CALL_MACS2_PEAKS.out.broad_igg
        peaks_broad_noigg = CALL_MACS2_PEAKS.out.broad_noigg
        peaks_broad_filtered = CALL_MACS2_PEAKS.out.broad_filtered
        peaks_broad_comb_igg_filtered = CALL_MACS2_PEAKS.out.broad_comb_igg_filtered
        peaks_broad_comb_igg = CALL_MACS2_PEAKS.out.broad_comb_igg
        ch_software_versions = ch_software_versions.mix(CALL_MACS2_PEAKS.out.versions)

    }


    /*
    * SUBWORKFLOW: QC peaks
    */
    ch_orig_csv = Channel.empty()
    ch_orig_widths = Channel.empty()
    ch_rip = Channel.empty()
    ch_rep_csv = Channel.empty()
    ch_conp_csv = Channel.empty()
    ch_conp_bed = Channel.empty()
    ch_conp_ann = Channel.empty()
    ch_orig_csv_comb_igg = Channel.empty()
    ch_orig_widths_comb_igg = Channel.empty()
    ch_rip_comb_igg = Channel.empty()
    ch_rep_csv_comb_igg = Channel.empty()
    ch_conp_csv_comb_igg = Channel.empty()
    ch_conp_bed_comb_igg = Channel.empty()
    ch_conp_ann_comb_igg = Channel.empty()
    if (params.run_peak_qc){
        QC_PEAKS(
            peaks_seacr_filtered,
            peaks_seacr_igg,
            peaks_seacr_noigg,
            peaks_seacr_comb_igg_filtered,
            peaks_seacr_comb_igg,
            peaks_narrow_filtered,
            peaks_narrow_igg,
            peaks_narrow_noigg,
            peaks_narrow_comb_igg_filtered,
            peaks_narrow_comb_igg,
            peaks_broad_filtered,
            peaks_broad_igg,
            peaks_broad_noigg,
            peaks_broad_comb_igg_filtered,
            peaks_broad_comb_igg,
            samplesheet,
            samplesheet_combine_igg,
            ch_samtools_bam,
            use_igg,
            params.min_replicates,
            params.fasta,
            params.gtf ? file(params.gtf, checkIfExists: true) : ch_dummy_file,
            srcdir
        )
        ch_orig_csv = QC_PEAKS.out.orig_csv
        ch_orig_widths = QC_PEAKS.out.orig_widths
        ch_rip = QC_PEAKS.out.rip
        ch_rep_csv = QC_PEAKS.out.rep_csv
        ch_conp_csv = QC_PEAKS.out.conp_csv
        ch_conp_bed = QC_PEAKS.out.conp_bed
        ch_conp_ann = QC_PEAKS.out.conp_ann
        ch_orig_csv_comb_igg = QC_PEAKS.out.orig_csv_comb_igg
        ch_orig_widths_comb_igg = QC_PEAKS.out.orig_widths_comb_igg
        ch_rip_comb_igg = QC_PEAKS.out.rip_comb_igg
        ch_rep_csv_comb_igg = QC_PEAKS.out.rep_csv_comb_igg
        ch_conp_csv_comb_igg = QC_PEAKS.out.conp_csv_comb_igg
        ch_conp_bed_comb_igg = QC_PEAKS.out.conp_bed_comb_igg
        ch_conp_ann_comb_igg = QC_PEAKS.out.conp_ann_comb_igg
        ch_software_versions = ch_software_versions.mix(QC_PEAKS.out.versions)
    }


    /*
    * SUBWORKFLOW: Call differential peaks from peaks using matched IgG controls
    */
    ch_dp = Channel.empty()
    ch_dp_comb_igg = Channel.empty()
    if (params.run_differential_peaks && params.comparison){
        CALL_DIFFERENTIAL_PEAKS(
            samplesheet,
            ch_conp_bed,
            samplesheet_combine_igg,
            ch_conp_bed_comb_igg,
            ch_samtools_bam,
            use_igg,
            srcdir
        )
        ch_dp = CALL_DIFFERENTIAL_PEAKS.out.dp
        ch_dp_comb_igg = CALL_DIFFERENTIAL_PEAKS.out.dp_comb_igg
        ch_software_versions = ch_software_versions.mix(CALL_DIFFERENTIAL_PEAKS.out.versions)
    }

    /*
    * SUBWORKFLOW: Generate analysis report
    */
    if (params.run_reporting){
        GENERATE_REPORT(
            ch_read_metrics,
            ch_frag_lens,
            samplesheet,
            ch_orig_csv,
            ch_orig_widths,
            ch_rip,
            ch_rep_csv,
            ch_conp_csv,
            ch_conp_bed,
            ch_conp_ann,
            ch_dp,
            samplesheet_combine_igg,
            ch_orig_csv_comb_igg,
            ch_orig_widths_comb_igg,
            ch_rip_comb_igg,
            ch_rep_csv_comb_igg,
            ch_conp_csv_comb_igg,
            ch_conp_bed_comb_igg,
            ch_conp_ann_comb_igg,
            ch_dp_comb_igg,
            use_igg,
            params.report_dir ? Channel.fromPath("${params.report_dir}", type: 'dir', checkIfExists: true) : Channel.fromPath("$projectDir/assets/local/report/", type: 'dir', checkIfExists: true),
            srcdir
        )
        ch_software_versions = ch_software_versions.mix(GENERATE_REPORT.out.versions)

    }


    /*
    * SUBWORKFLOW: make heatmaps and other summary plots for peaks using matched IgG controls
    */
    if (params.run_summary_plots){
        SUMMARY_PLOTS(
            ch_bigwig_markdup,
            ch_conp_bed,
            ch_conp_ann,
            ch_dp,
            ch_conp_bed_comb_igg,
            ch_conp_ann_comb_igg,
            ch_dp_comb_igg,
            use_igg,
            params.gtf ? file(params.gtf, checkIfExists: true) : ch_dummy_file,
            params.gene_bed ? file(params.gene_bed, checkIfExists: true) : "$projectDir/assets/dummy_file.txt"
        )
        ch_software_versions = ch_software_versions.mix(SUMMARY_PLOTS.out.versions)
    }


    /*
    * "process control" workflow: peak calling, QC and analysis report
    */
    ch_peaks_all = Channel.empty()
    ch_peaks_final = Channel.empty()
    if (params.workflow == "process_controls"){
        /*
        * SUBWORKFLOW: Call peaks for "process controls" workflow
        */
        if(params.run_peak_calling) {
            CALL_PEAKS_PROCESS_CONTROLS (
            ch_bedgraph_markdup,
            ch_bedgraph_dedup,
            ch_samtools_bam_markdup,
            use_igg,
            params.igg_dir
            )
            ch_peaks_all = CALL_PEAKS_PROCESS_CONTROLS.out.peaks_all
            // ch_peaks_all.view()
            // [ meta, [peaks] ]

            ch_peaks_final = CALL_PEAKS_PROCESS_CONTROLS.out.peaks_final
            // ch_peaks_final.view()
            // [ meta, [peaks] ]

            ch_software_versions = ch_software_versions.mix(CALL_PEAKS_PROCESS_CONTROLS.out.versions)
        }

        /*
        * SUBWORKFLOW: QC peaks for "process controls" workflow
        */
        if (params.run_peak_qc){
            ch_ref_peaks = params.ref_peaks ? Channel.fromPath("${params.ref_peaks}", type: "dir", checkIfExists: true) : Channel.empty()
            QC_PROCESS_CONTROLS(
            samplesheet,
            params.genome,
            ch_samtools_bam,
            ch_peaks_all,
            ch_peaks_final,
            ch_ref_peaks
            )
            ch_software_versions = ch_software_versions.mix(QC_PROCESS_CONTROLS.out.versions)
        }

        /*
        * SUBWORKFLOW: Generate report for "process controls" workflow
        */
        if (params.run_reporting){
            GENERATE_REPORT_PROCESS_CONTROLS(
                samplesheet,
                ch_read_metrics,
                ch_frag_lens.collect{it[1]},
                QC_PROCESS_CONTROLS.out.orig_csv.collect{it[1]},
                QC_PROCESS_CONTROLS.out.orig_widths.collect{it[1]},
                QC_PROCESS_CONTROLS.out.rip.collect{it[1]},
                QC_PROCESS_CONTROLS.out.rep_csv.collect{it[1]},
                params.saved_data ? Channel.fromPath("${params.saved_data}", type: "dir", checkIfExists: true) : Channel.empty(),
                params.report_rmd ? Channel.fromPath("${params.report_rmd}", checkIfExists: true) : Channel.empty()
            )
            ch_software_versions = ch_software_versions.mix(GENERATE_REPORT_PROCESS_CONTROLS.out.versions)

        }

    }

    /*
    * Write output csv files with paths to important intermediate files for downstream analysis
    */
    if (params.update == 'output'){
        /*
        UPDATE_OUTPUT(
            samplesheet,
            srcdir,
            outdir
        )
        */
    }else{

        WRITE_OUTPUT_CSV(
            ch_samtools_bam_markdup.join(ch_samtools_bai_markdup),
            ch_bedgraph_markdup,
            ch_bedgraph_dedup,
            peaks_narrow_noigg,
            peaks_narrow_filtered,
            peaks_narrow_igg,
            peaks_narrow_comb_igg_filtered,
            peaks_narrow_comb_igg,
            peaks_broad_noigg,
            peaks_broad_filtered,
            peaks_broad_igg,
            peaks_broad_comb_igg_filtered,
            peaks_broad_comb_igg,
            peaks_seacr_noigg,
            peaks_seacr_filtered,
            peaks_seacr_igg,
            peaks_seacr_comb_igg_filtered,
            peaks_seacr_comb_igg,
            ch_conp_bed,
            ch_conp_bed_comb_igg
        )


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
