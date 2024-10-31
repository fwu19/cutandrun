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
    params.input
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check spike-in
checkPathParamList = [
        params.spikein_bowtie2,
        params.spikein_fasta
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Check mandatory parameters that cannot be checked in the groovy lib as we want a channel for them
if (params.input) { ch_input = file(params.input) } else { exit 1, "Input samplesheet not specified!" }

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
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

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
include { INPUT_CHECK                } from "../subworkflows/local/input_check"
include { CUT as PEAK_TO_BED         } from '../modules/local/linux/cut'
include { AWK as AWK_NAME_PEAK_BED   } from "../modules/local/linux/awk"
include { IGV_SESSION                } from "../modules/local/python/igv_session"
include { AWK as AWK_EXTRACT_SUMMITS } from "../modules/local/linux/awk"
include { SAMTOOLS_CUSTOMVIEW        } from "../modules/local/samtools_custom_view"
include { FRAG_LEN_HIST              } from "../modules/local/python/frag_len_hist"
//include { MULTIQC                    } from "../modules/local/multiqc"
include { BEDTOOLS_INTERSECT as SEACR_PEAKS_BEDTOOLS_INTERSECT   } from "../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT   } from "../modules/nf-core/bedtools/intersect/main"
include { BEDTOOLS_INTERSECT as MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT   } from "../modules/nf-core/bedtools/intersect/main"

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
include { CONSENSUS_PEAKS                                  } from "../subworkflows/local/consensus_peaks"
include { CONSENSUS_PEAKS as CONSENSUS_PEAKS_ALL           } from "../subworkflows/local/consensus_peaks"
include { EXTRACT_FRAGMENTS                                } from "../subworkflows/local/extract_fragments"
include { COMPUTE_GENOMECOVERAGE                           } from "../subworkflows/local/compute_genomecoverage"
include { DEEPTOOLS_QC                                     } from "../subworkflows/local/deeptools_qc"
include { PEAK_QC                                          } from "../subworkflows/local/peak_qc"
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
include { SEACR_CALLPEAK as SEACR_CALLPEAK_IGG                         } from "../modules/nf-core/seacr/callpeak/main"
include { SEACR_CALLPEAK as SEACR_CALLPEAK_NOIGG                       } from "../modules/nf-core/seacr/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_IGG_NARROW                         } from "../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NOIGG_NARROW                       } from "../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_IGG_BROAD                         } from "../modules/nf-core/macs2/callpeak/main"
include { MACS2_CALLPEAK as MACS2_CALLPEAK_NOIGG_BROAD                       } from "../modules/nf-core/macs2/callpeak/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE      } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS     } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE          } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS         } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_GENE_ALL  } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_COMPUTEMATRIX as DEEPTOOLS_COMPUTEMATRIX_PEAKS_ALL } from "../modules/nf-core/deeptools/computematrix/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_GENE_ALL      } from "../modules/nf-core/deeptools/plotheatmap/main"
include { DEEPTOOLS_PLOTHEATMAP as DEEPTOOLS_PLOTHEATMAP_PEAKS_ALL     } from "../modules/nf-core/deeptools/plotheatmap/main"
include { CUSTOM_DUMPSOFTWAREVERSIONS                                  } from "../modules/local/custom_dumpsoftwareversions"

/*
 * SUBWORKFLOWS
 */



/*
========================================================================================
    IMPORT CUSTOM MODULES/SUBWORKFLOWS
========================================================================================
*/

include { MULTIQC  } from '../modules/local2/multiqc'
include { FRAGMENT_LENGTH  } from '../modules/local2/fragment_length'
include { READS_IN_PEAK  } from '../modules/local2/reads_in_peak'
include { READ_METRICS  } from '../modules/local2/read_metrics'
include { ORIGINAL_PEAK_METRICS  } from '../modules/local2/original_peak_metrics'
include { REPLICATED_PEAKS } from '../modules/local2/replicated_peaks'
include { CONSENSUS_PEAKS } from '../modules/local2/consensus_peaks'
include { PLOT_METRICS } from '../modules/local2/plot_metrics'



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
        INPUT_CHECK (
            ch_input
        )

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

    ch_bedgraph               = Channel.empty()
    ch_bigwig                 = Channel.empty()
    ch_seacr_peaks_igg            = Channel.empty()
    ch_seacr_peaks_noigg            = Channel.empty()
    ch_seacr_peaks            = Channel.empty()
    ch_macs2_peaks            = Channel.empty()
    ch_macs2_peaks_igg_narrow = Channel.empty()
    ch_macs2_peaks_noigg_narrow = Channel.empty()
    ch_macs2_peaks_narrow     = Channel.empty()
    ch_macs2_peaks_broad      = Channel.empty()
    ch_macs2_peaks_igg_broad      = Channel.empty()
    ch_macs2_peaks_noigg_broad      = Channel.empty()
    ch_peaks_primary          = Channel.empty()
    ch_peaks_secondary        = Channel.empty()
    ch_peaks_summits          = Channel.empty()
    ch_consensus_peaks        = Channel.empty()
    ch_consensus_peaks_unfilt = Channel.empty()
    if(params.run_alignment && params.run_read_filter) {
        /*
        * SUBWORKFLOW: Convert BAM files to bedgraph/bigwig and apply spikein normalisation if required
        */
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

    if(params.run_peak_calling) {

        /*
         * CHANNEL: Separate bedgraphs into target/control for SEACR
         * for SEACR, use markdup for target and dedup for control
         */
        ch_bedgraph_markdup
            .filter { it -> it[0].is_control == false }
            .set { ch_bedgraph_target }
        ch_bedgraph_dedup
            .filter { it -> it[0].is_control == true }
            .set { ch_bedgraph_control }
        //ch_bedgraph_target | view
        //ch_bedgraph_control | view


        /*
        * CHANNEL: Separate bams into target/control
        * for MACS2 use markdup for both target and control
        */
        ch_samtools_bam_markdup
            .filter { it -> it[0].is_control == false }
            .set { ch_bam_target }
        ch_samtools_bam_markdup
            .filter { it -> it[0].is_control == true }
            .set { ch_bam_control }
        //ch_bam_target | view
        //ch_bam_control | view

        if(params.use_control) {
            /*
            * MODULE: Call peaks using SEACR with IgG control
            */
            if('seacr' in callers) {
                /*
                * CHANNEL: Create target/control pairings
                */
                ch_bedgraph_control.map{ row -> [row[0].control_group, row] }
                .cross( ch_bedgraph_target.map{ row -> [row[0].control_group, row] } )
                .map {
                    row ->
                    [ row[1][1][0], row[1][1][1], row[0][1][1] ]
                }
                .set { ch_bedgraph_paired }
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BEDGRAPH, CONTROL_BEDGRAPH]

                SEACR_CALLPEAK_IGG (
                    ch_bedgraph_paired,
                    params.seacr_peak_threshold
                )
                ch_seacr_peaks_igg       = SEACR_CALLPEAK_IGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_IGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_CALLPEAK_IGG.out.bed | view
            }

            if('seacr' in callers) {
                /*
                * CHANNEL: Add fake control channel
                */
                ch_bedgraph_target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_bedgraph_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
                // ch_bedgraph_target_fctrl | view

                SEACR_CALLPEAK_NOIGG (
                    ch_bedgraph_target_fctrl,
                    params.seacr_peak_threshold
                )
                ch_seacr_peaks_noigg       = SEACR_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_NO_IGG.out.bed | view
            }

            if('seacr' in callers) {
                /*
                * CHANNEL: mix igg and noigg SEACR peaks
                */

                ch_seacr_peaks_intersect = ch_seacr_peaks_igg
                    .join (ch_seacr_peaks_noigg)
                    .map {row -> [row[0], row[1], row[2]]}

                SEACR_PEAKS_BEDTOOLS_INTERSECT(
                    ch_seacr_peaks_intersect,
                    [[:],[]]
                )
                ch_seacr_peaks = SEACR_PEAKS_BEDTOOLS_INTERSECT.out.intersect
                ch_software_versions = ch_software_versions.mix(SEACR_PEAKS_BEDTOOLS_INTERSECT.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_PEAKS_BEDTOOLS_INTERSECT.out.intersect | view
            }

            if('macs2' in callers) {
                /*
                * CHANNEL: Create target/control pairings
                */
                ch_bam_control.map{ row -> [row[0].control_group, row] }
                .cross( ch_bam_target.map{ row -> [row[0].control_group, row] } )
                .map {
                    row ->
                    [ row[1][1][0], row[1][1][1], row[0][1][1] ]
                }
                .set { ch_bam_paired }
                // EXAMPLE CHANNEL STRUCT: [[META], TARGET_BAM, CONTROL_BAM]
                //ch_bam_paired | view

                MACS2_CALLPEAK_IGG_NARROW (
                    ch_bam_paired,
                    params.macs_gsize
                )
                ch_macs2_peaks_igg_narrow       = MACS2_CALLPEAK_IGG_NARROW.out.peak
                ch_peaks_summits_igg_narrow     = MACS2_CALLPEAK_IGG_NARROW.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_IGG_NARROW.out.versions)

                MACS2_CALLPEAK_IGG_BROAD (
                    ch_bam_paired,
                    params.macs_gsize
                )
                ch_macs2_peaks_igg_broad       = MACS2_CALLPEAK_IGG_BROAD.out.peak
                ch_peaks_summits_igg_broad     = MACS2_CALLPEAK_IGG_BROAD.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_IGG_BROAD.out.versions)

                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //MACS2_CALLPEAK_IGG.out.peak | view
            }

            if('macs2' in callers) {

                /*
                * CHANNEL: Add fake control channel
                */
                ch_bam_target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_samtools_bam_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]
                //ch_samtools_bam_target_fctrl | view

                MACS2_CALLPEAK_NOIGG_NARROW (
                    ch_samtools_bam_target_fctrl,
                    params.macs_gsize
                )
                ch_macs2_peaks_noigg_narrow       = MACS2_CALLPEAK_NOIGG_NARROW.out.peak
                ch_peaks_summits_noigg_narrow     = MACS2_CALLPEAK_NOIGG_NARROW.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_NOIGG_NARROW.out.versions)

                MACS2_CALLPEAK_NOIGG_BROAD (
                    ch_samtools_bam_target_fctrl,
                    params.macs_gsize
                )
                ch_macs2_peaks_noigg_broad       = MACS2_CALLPEAK_NOIGG_BROAD.out.peak
                ch_peaks_summits_noigg_broad     = MACS2_CALLPEAK_NOIGG_BROAD.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_NOIGG_NARROW.out.versions)

                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                // MACS2_CALLPEAK_NOIGG.out.peak | view
            }

            if('macs2' in callers) {
                /*
                * CHANNEL: mix igg and noigg MACS2 narrow peaks
                */

                ch_macs2_peaks_narrow_intersect = ch_macs2_peaks_igg_narrow
                    .join (ch_macs2_peaks_noigg_narrow)
                    .map {row -> [row[0], row[1], row[2]]}

                MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT(
                    ch_macs2_peaks_narrow_intersect,
                    [[:],[]]
                )
                ch_macs2_peaks_narrow = MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.intersect
                ch_software_versions = ch_software_versions.mix(MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //MACS2_PEAKS_NARROW_BEDTOOLS_INTERSECT.out.intersect | view

                /*
                * CHANNEL: mix igg and noigg MACS2 broad peaks
                */

                ch_macs2_peaks_broad_intersect = ch_macs2_peaks_noigg_broad
                    .join (ch_macs2_peaks_igg_broad)
                    .map {row -> [row[0], row[1], row[2]]}

                MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT(
                    ch_macs2_peaks_broad_intersect,
                    [[:],[]]
                )
                ch_macs2_peaks_broad = MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.intersect
                ch_software_versions = ch_software_versions.mix(MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //MACS2_PEAKS_BROAD_BEDTOOLS_INTERSECT.out.intersect | view

            }

        }
        else {
            /*
            * MODULE: Call peaks without IgG Control
            */
            if('seacr' in callers) {
                /*
                * CHANNEL: Add fake control channel
                */
                ch_bedgraph_target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_bedgraph_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BED, FAKE_CTRL]
                // ch_bedgraph_target_fctrl | view

                SEACR_CALLPEAK_NOIGG (
                    ch_bedgraph_target_fctrl,
                    params.seacr_peak_threshold
                )
                ch_seacr_peaks       = SEACR_CALLPEAK_NOIGG.out.bed
                ch_software_versions = ch_software_versions.mix(SEACR_CALLPEAK_NOIGG.out.versions)
                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                //SEACR_NO_IGG.out.bed | view
            }

            if('macs2' in callers) {
                /*
                * CHANNEL: Add fake control channel
                */
                ch_bam_target.map{ row-> [ row[0], row[1], [] ] }
                .set { ch_samtools_bam_target_fctrl }
                // EXAMPLE CHANNEL STRUCT: [[META], BAM, FAKE_CTRL]
                //ch_samtools_bam_target_fctrl | view

                MACS2_CALLPEAK_NOIGG_NARROW (
                    ch_samtools_bam_target_fctrl,
                    params.macs_gsize
                )
                ch_macs2_peaks_narrow       = MACS2_CALLPEAK_NOIGG_NARROW.out.peak
                ch_peaks_summits_narrow     = MACS2_CALLPEAK_NOIGG_NARROW.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_NOIGG_NARROW.out.versions)

                MACS2_CALLPEAK_NOIGG_BROAD (
                    ch_samtools_bam_target_fctrl,
                    params.macs_gsize
                )
                ch_macs2_peaks_broad       = MACS2_CALLPEAK_NOIGG_BROAD.out.peak
                ch_peaks_summits_broad     = MACS2_CALLPEAK_NOIGG_BROAD.out.bed
                ch_software_versions = ch_software_versions.mix(MACS2_CALLPEAK_NOIGG_BROAD.out.versions)

                // EXAMPLE CHANNEL STRUCT: [[META], BED]
                // MACS2_CALLPEAK_NOIGG.out.peak | view
            }
        }

    }


    if (params.run_local){
        /* Run MultiQC */
        MULTIQC(
            ch_multiqc_custom_config,
            ch_bowtie2_log.collect{it[1]}.ifEmpty([]),
            ch_bowtie2_spikein_log.collect{it[1]}.ifEmpty([]),
            ch_samtools_stats.collect{it[1]}.ifEmpty([]),
            ch_samtools_flagstat.collect{it[1]}.ifEmpty([]),
            ch_samtools_idxstats.collect{it[1]}.ifEmpty([]),
            ch_markduplicates_metrics.collect{it[1]}.ifEmpty([])

        )

        /* Compute fragment length */
        FRAGMENT_LENGTH(ch_samtools_bam)
        ch_frag_len = FRAGMENT_LENGTH.out.frag_len
            .flatten()
            .collect()
        //ch_frag_len.view()

        /* Compute reads in peak */
        ch_macs2_peaks_narrow
            .concat(
                ch_macs2_peaks_igg_narrow,
                ch_macs2_peaks_noigg_narrow,
                ch_macs2_peaks_broad,
                ch_macs2_peaks_igg_broad,
                ch_macs2_peaks_noigg_broad,
                ch_seacr_peaks,
                ch_seacr_peaks_igg,
                ch_seacr_peaks_noigg
            )
            .groupTuple(by: 0)
            .join(ch_samtools_bam)
            .set{ch_peak_bam}

        READS_IN_PEAK( ch_peak_bam )
        rip_list = READS_IN_PEAK.out.rip
            .map{ it -> it[1] }
            .flatten()
            .collect()
        //ch_rip.view()

        /* Collect reads metrics */
        READ_METRICS( params.input, MULTIQC.out.data, ch_frag_len )
        read_metrics = READ_METRICS.out.metrics
        //read_metrics.view()

        /* Collect metrics for original peaks */
        ch_macs2_peaks_narrow
            .concat(
                ch_macs2_peaks_igg_narrow,
                ch_macs2_peaks_noigg_narrow,
                ch_macs2_peaks_broad,
                ch_macs2_peaks_igg_broad,
                ch_macs2_peaks_noigg_broad,
                ch_seacr_peaks,
                ch_seacr_peaks_igg,
                ch_seacr_peaks_noigg
            )
            .map{ it -> it[1]  }
            .collect()
            .set{ peak_list }

        ORIGINAL_PEAK_METRICS( read_metrics, peak_list, rip_list)
        original_peaks = ORIGINAL_PEAK_METRICS.out.metrics
        // original_peaks.view()

        /* Generate replicated peaks and collect metrics  */
        REPLICATED_PEAKS( params.input, original_peaks, params.min_replicates )
        replicated_peaks = REPLICATED_PEAKS.out.metrics

        /* Generate consensus peaks and collect metrics  */
        CONSENSUS_PEAKS( replicated_peaks )
        consensus_peaks = CONSENSUS_PEAKS.out.metrics

        /* Make plots for report */
        PLOT_METRICS( read_metrics, original_peaks, replicated_peaks, consensus_peaks )
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
