process READS_IN_PEAK {
    label "process_single"
    tag "READS_IN_PEAK on ${meta.id}"

//    conda "bioconda::samtools=1.17 bioconda::bedtools=2.30.0"
    module = [ 'SAMtools/1.11-GCC-10.2.0', 'BEDTools/2.30.0-GCC-10.2.0' ]

    input:
    tuple val( meta ), path( "peaks/*" ), path(bam)

    output:
    tuple val( meta ), path( "*.reads_in_peak.csv" ), emit: rip

    script:
    """
    reads_in_peak.sh ${meta.id} $bam
    """
}
