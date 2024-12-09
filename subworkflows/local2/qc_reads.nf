/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { FRAGMENT_LENGTHS           } from '../../modules/local2/fragment_lengths'
include { READ_METRICS              } from '../../modules/local2/read_metrics'

workflow QC_READS {
    take:
    ch_multiqc_data
    ch_samtools_bam


    main:
    ch_read_metrics = Channel.empty()
    ch_frag_lens = Channel.empty()
    /*
    * Collect reads metrics from MultiQC data
    */

    READ_METRICS(
        ch_multiqc_data
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


    emit:

    read_metrics = ch_read_metrics
    frag_lens = ch_frag_lens

}
