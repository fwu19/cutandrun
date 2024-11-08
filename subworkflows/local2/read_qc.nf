/*
* Collect read and peak metrics and generate replicated and consensus peaks
*/

include { FRAGMENT_LENGTH           } from '../../modules/local2/fragment_length'
include { READ_METRICS              } from '../../modules/local2/read_metrics'

workflow READ_QC {
    take:
    multiqc
    bam


    main:

    /* Compute fragment length */
    FRAGMENT_LENGTH(
        bam
    )
    //ch_frag_len.view()


    /* Collect reads metrics */
    READ_METRICS(
        multiqc,
        FRAGMENT_LENGTH.out.frag_len.flatten().collect()
    )
    //READ_METRICS.out.metrics.view()

    emit:
    metrics = READ_METRICS.out.data // *.csv

}
