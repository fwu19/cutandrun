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
    params.blacklist,
    params.bowtie2,
    params.fasta,
    params.gtf,
    //params.input
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
ch_dummy_csv = file("$projectDir/assets/local/dummy_file.csv", checkIfExists: true)

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

// Header files for MultiQC
ch_frag_len_header_multiqc              = file("$projectDir/assets/multiqc/frag_len_header.txt", checkIfExists: true)
ch_frip_score_header_multiqc            = file("$projectDir/assets/multiqc/frip_score_header.txt", checkIfExists: true)
ch_peak_counts_header_multiqc           = file("$projectDir/assets/multiqc/peak_counts_header.txt", checkIfExists: true)
ch_peak_counts_consensus_header_multiqc = file("$projectDir/assets/multiqc/peak_counts_consensus_header.txt", checkIfExists: true)
ch_peak_reprod_header_multiqc           = file("$projectDir/assets/multiqc/peak_reprod_header.txt", checkIfExists: true)
ch_linear_duplication_header_multiqc    = file("$projectDir/assets/multiqc/linear_duplication_header.txt", checkIfExists: true)


/*
========================================================================================
    INIALISE PARAMETERS AND OPTIONS
========================================================================================
*/

// Init aligners
def prepare_tool_indices = ["bowtie2"]

// Check peak caller params
def caller_list = ['seacr', 'macs2']
callers = params.peakcaller ? params.peakcaller.split(',').collect{ it.trim().toLowerCase() } : ['seacr']
if ((caller_list + callers).unique().size() != caller_list.size()) {
    exit 1, "Invalid variant calller option: ${params.peakcaller}. Valid options: ${caller_list.join(', ')}"
}

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { CUT as PEAK_TO_BED         } from '../modules/local/linux/cut'
include { AWK as AWK_NAME_PEAK_BED   } from "../modules/local/linux/awk"
include { IGV_SESSION                } from "../modules/local/python/igv_session"
include { AWK as AWK_EXTRACT_SUMMITS } from "../modules/local/linux/awk"
include { SAMTOOLS_CUSTOMVIEW        } from "../modules/local/samtools_custom_view"
include { FRAG_LEN_HIST              } from "../modules/local/python/frag_len_hist"
//include { MULTIQC                    } from "../modules/local/multiqc"

/*
 * SUBWORKFLOWS
 */
include { PREPARE_GENOME                                   } from "../subworkflows/local/prepare_genome"
include { FASTQC_TRIMGALORE                                } from "../subworkflows/local/fastqc_trimgalore"
include { ALIGN_BOWTIE2                                    } from "../subworkflows/local/align_bowtie2"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_TARGET_META  } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_BT2_SPIKEIN_META } from "../subworkflows/local/extract_metadata_awk"
include { EXTRACT_METADATA_AWK as EXTRACT_PICARD_DUP_META  } from "../subworkflows/local/extract_metadata_awk"
include { MARK_DUPLICATES_PICARD                           } from "../subworkflows/local/mark_duplicates_picard"
include { MARK_DUPLICATES_PICARD as DEDUPLICATE_PICARD     } from "../subworkflows/local/mark_duplicates_picard"
include { EXTRACT_FRAGMENTS                                } from "../subworkflows/local/extract_fragments"
include { SAMTOOLS_VIEW_SORT_STATS as FILTER_READS         } from "../subworkflows/local/samtools_view_sort_stats"
include { DEDUPLICATE_LINEAR                               } from "../subworkflows/local/deduplicate_linear"

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

/*
 * MODULES
 */
include { CAT_FASTQ                                                    } from "../modules/nf-core/cat/fastq/main"
include { PRESEQ_LCEXTRAP                                              } from "../modules/nf-core/preseq/lcextrap/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS                                  } from "../modules/local/custom_dumpsoftwareversions"

/*
 * SUBWORKFLOWS
 */



/*
========================================================================================
    IMPORT CUSTOM MODULES/SUBWORKFLOWS
========================================================================================
*/
include { INPUT_CHECK                   } from "../subworkflows/local2/input_check"
include { CALL_PEAKS                    } from '../subworkflows/local2/call_peaks'
include { COMPUTE_GENOMECOVERAGE        } from "../subworkflows/local2/compute_genomecoverage"

include { GET_FASTQ_PATHS               } from '../modules/local2/get_fastq_paths'
include { MULTIQC                       } from '../modules/local2/multiqc'
include { FRAGMENT_LENGTHS               } from '../modules/local2/fragment_lengths'
include { READ_METRICS                  } from '../modules/local2/read_metrics'
include { READS_IN_PEAK                 } from '../modules/local2/reads_in_peak'
include { ORIGINAL_PEAKS                } from '../modules/local2/original_peaks'
include { ORIGINAL_PEAK_WIDTHS          } from '../modules/local2/original_peak_widths'
include { REPLICATED_PEAKS              } from '../modules/local2/replicated_peaks'
include { CONSENSUS_PEAKS               } from '../modules/local2/consensus_peaks'
include { ANNOTATE_CONSENSUS_PEAKS      } from '../modules/local2/annotate_consensus_peaks'
include { GENERATE_REPORT                  } from '../modules/local2/generate_report'
include { READS_IN_CONSENSUS_PEAKS      } from '../modules/local2/reads_in_consensus_peaks'
include { DIFFERENTIAL_PEAKS            } from '../modules/local2/differential_peaks'



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
        
        if ( params.input_dir =~ 'dummy' ){
            if ( params.input =~ 'dummy' ){
                exit 1, 'Neither --input nor --input_dir is specified!'
            }else {
                ch_input = Channel.fromPath( params.input, checkIfExists: true )
            }
        }else {
            GET_FASTQ_PATHS (
                params.input_dir
            )
            ch_input = GET_FASTQ_PATHS.out.csv
        }

        ch_metadata = params.metadata ? file( params.metadata, checkIfExists: true ) : ch_dummy_csv
        INPUT_CHECK (
            ch_input,
            ch_metadata
        )

        samplesheet = INPUT_CHECK.out.samplesheet

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
            ch_software_versions          = ch_software_versions.mix(ALIGN_BOWTIE2.out.versions)
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
     * MODULE: Run preseq on BAM files before de-duplication
    */
    ch_preseq_output = Channel.empty()
    if (params.run_preseq) {
        PRESEQ_LCEXTRAP (
            ch_samtools_bam
        )
        ch_preseq_output = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_software_versions = ch_software_versions.mix(PRESEQ_LCEXTRAP.out.versions)
    }

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
     * SUBWORKFLOW: Remove linear amplification duplicates - default is false
     */
    ch_linear_metrics         = Channel.empty()
    ch_linear_duplication_mqc = Channel.empty()
    if (params.run_remove_linear_dups) {
        DEDUPLICATE_LINEAR (
            ch_samtools_bam,
            ch_samtools_bai,
            PREPARE_GENOME.out.fasta.collect(),
            PREPARE_GENOME.out.fasta_index.collect(),
            params.dedup_target_reads,
            ch_linear_duplication_header_multiqc
        )
        ch_samtools_bam           = DEDUPLICATE_LINEAR.out.bam
        ch_samtools_bai           = DEDUPLICATE_LINEAR.out.bai
        ch_samtools_stats         = DEDUPLICATE_LINEAR.out.stats
        ch_samtools_flagstat      = DEDUPLICATE_LINEAR.out.flagstat
        ch_samtools_idxstats      = DEDUPLICATE_LINEAR.out.idxstats
        ch_linear_metrics         = DEDUPLICATE_LINEAR.out.metrics
        ch_linear_duplication_mqc = DEDUPLICATE_LINEAR.out.linear_metrics_mqc
        ch_software_versions      = ch_software_versions.mix(DEDUPLICATE_LINEAR.out.versions)
    }


    /*
    * SUBWORKFLOW: Convert BAM files to bedgraph/bigwig and apply spikein normalisation if required
    */
    ch_bedgraph_markdup     = Channel.empty()
    ch_bedgraph_dedup       = Channel.empty()
    if(params.run_alignment && params.run_read_filter) {
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
        //ch_bigwig_markdup            = COMPUTE_GENOMECOVERAGE_MARKDUP.out.bigwig_cpm
        ch_software_versions = ch_software_versions.mix(COMPUTE_GENOMECOVERAGE.out.versions)

    }

    /*
    * Run MultiQC
    */
    multiqc_data = Channel.empty()
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

    }

    /*
     * SUBWORKFLOW: Call peaks from individual samples
     */
    ch_peaks_all = Channel.empty()
    ch_peaks_final = Channel.empty()
    if(params.run_peak_calling) {
        CALL_PEAKS (
            ch_bedgraph_markdup,
            ch_bedgraph_dedup,
            ch_samtools_bam_markdup
        )
        ch_peaks_all = CALL_PEAKS.out.peaks_all
        // ch_peaks_all.view()
        // [ meta, [peaks] ]

        ch_peaks_final = CALL_PEAKS.out.peaks_final
        // ch_peaks_final.view()
        // [ meta, [peaks] ]
    }

    ch_read_metrics = Channel.empty()
    ch_frag_lens = Channel.empty()
    if (params.run_local_read_qc){

        /*
        * Collect reads metrics from MultiQC data
        */

        READ_METRICS(
            MULTIQC.out.data
        )
        ch_read_metrics = READ_METRICS.out.csv
        //ch_read_metrics.view()
        // path(csv)

        /*
        * Compute fragment length
        */

        FRAGMENT_LENGTHS(
            ch_samtools_bam
        )
        ch_frag_lens = FRAGMENT_LENGTHS.out.txt
        // ch_frag_lens.view()
        // [ meta, path(txt) ]
    }

    ch_rip = Channel.empty()
    ch_orig_csv = Channel.empty()
    ch_orig_widths = Channel.empty()
    ch_rep_bed = Channel.empty()
    ch_rep_csv = Channel.empty()
    ch_con_bed = Channel.empty()
    ch_con_csv = Channel.empty()
    ch_conp_ann = Channel.empty()
    if (params.run_local_peak_qc){
        /*
        * Compute reads in peak
        */
        ch_peaks_final
            .join(ch_samtools_bam)
            .set{ch_peak_bam}

        READS_IN_PEAK(
            ch_peak_bam
        )
        ch_rip = READS_IN_PEAK.out.csv.collect{it[1]}
        // ch_rip.view()


        /*
        * Collect metrics for original peaks
        */
        ORIGINAL_PEAKS(
            ch_peaks_all
        )
        ch_orig_csv = ORIGINAL_PEAKS.out.csv
        // ch_orig_peaks.view()
        // path(peak_metrics)

        /*
        * Collect peak widths for final peaks
        */

        ORIGINAL_PEAK_WIDTHS(
            ch_peaks_final
        )
        ch_orig_widths = ORIGINAL_PEAK_WIDTHS.out.csv
        // ch_orig_widths.view()
        // path(peak_widths)

        /*
        * Generate replicated peaks and collect metrics
        */
        REPLICATED_PEAKS(
            ch_peaks_final
            .filter { it[0].call_rep_peak == true }
            .map { it -> [ [it[0].group, it[0].target], it[1] ]}
            .groupTuple (by: 0)
            .map { it -> [ it[0][0], it[0][1], it[1].flatten().collect() ] },
            params.min_replicates
        )
        ch_rep_bed = REPLICATED_PEAKS.out.bed
        ch_rep_csv = REPLICATED_PEAKS.out.csv
        // ch_rep_bed.view()
        // [ target, [peaks] ]

        /*
        * Generate consensus peaks and collect metrics
        */
        CONSENSUS_PEAKS(
            samplesheet,
            ch_rep_bed
            .groupTuple ( by: 0 )
            .map { it -> [ it[0], it[1].flatten().collect() ] }
        )
        ch_con_bed = CONSENSUS_PEAKS.out.bed
        ch_con_csv = CONSENSUS_PEAKS.out.csv
        // ch_conp_bed.view()
        // [ target, path(conp) ]

        /*
        * Annotate consensus peaks
        */
        ch_gtf_ann = params.local_assets ? file("${params.local_assets}/${params.genome}/genes.proteinCoding_lncRNA.gtf") : ch_dummy_file

        ANNOTATE_CONSENSUS_PEAKS(
            params.genome,
            ch_gtf_ann,
            ch_con_bed.collect{it[1]}.flatten()
        )
        ch_conp_ann = ANNOTATE_CONSENSUS_PEAKS.out.txt
        // ch_conp_ann.view()

    }

    ch_conp_reads = Channel.empty()
    ch_dp = Channel.empty()
    if (params.run_local_dp){
            /*
            * Count reads in consensus peaks
            */
            // ch_samtools_bam.view()
            ch_con_peaks
                .cross (
                    ch_samtools_bam
                        .filter { it[0].target != "IgG" }
                        .map { it -> [ it[0].target, it ] }
                )
                .map { it -> [ it[1][1][0], it[1][1][1], it[0][1] ] }
                .set { ch_bam_conp }
            //ch_bam_conp.view()
            // [ meta, bam, [conp] ]


            READS_IN_CONSENSUS_PEAKS(
                ch_bam_conp
            )
            ch_conp_reads = READS_IN_CONSENSUS_PEAKS.out.count
            // ch_conp_reads.view()
            // [ target, [path/to/fragmentCounts.txt] ]

            ch_conp_reads
                .groupTuple( by: 0 )
                .map { it -> [ it[0], it[1].flatten().collect() ]}
                .cross ( CONSENSUS_PEAKS.out.bed )
                .map { it -> [ it[0][0], it[0][1], it[1][1] ]}
                .set { ch_tgt_reads_conp }
            //ch_tgt_reads_conp.view()
            // [ target, [read_count], [conp] ]

            /*
            * Call differential peaks
            */
            ch_comparison = params.comparison ? file ( params.comparison, checkIfExists: true ) : ch_dummy_file

            DIFFERENTIAL_PEAKS(
                samplesheet,
                ch_comparison,
                ch_tgt_reads_conp
            )
            ch_dp = DIFFERENTIAL_PEAKS.out.data
            //ch_dp.view()
            // [ path(*.{rds,csv}) ]

    }

    if (params.run_local_report){
            /*
            * Make plots for report
            */
            GENERATE_REPORT(
                samplesheet,
                ch_read_metrics.ifEmpty([]),
                ch_frag_lens.ifEmpty([]),
                ch_orig_csv.collect{it[1]}.ifEmpty([]),
                ch_orig_widths.collect{it[1]}.ifEmpty([]),
                ch_rip.collect{it[1]}.ifEmpty([]),
                ch_rep_csv.collect{it[1]}.ifEmpty([]),
                ch_con_csv.collect{it[1]}.ifEmpty([]),
                ch_con_bed.collect{it[1]}.flatten().collect().ifEmpty([]),
                ch_conp_ann.collect{it[1]}.ifEmpty([]),
                ch_dp.flatten().collect().ifEmpty([]),
                file("$projectDir/assets/local/report.Rmd", checkIfExists: true)
            )


    }



}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
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
