process READS_IN_PEAK {
    label "process_single"
    tag "READS_IN_PEAK on ${meta.id}"

//    conda "bioconda::samtools=1.17 bioconda::bedtools=2.30.0"
    module = [ 'SAMtools/1.11-GCC-10.2.0', 'BEDTools/2.30.0-GCC-10.2.0' ]

    input:
    tuple val( meta ), path( "peaks/*" ), path(bam)

    output:
    path( "*.reads_in_peak.csv" ), emit: csv
    path ( "versions.yml" ), emit: versions

    script:
    """
    reads_in_peak.sh ${meta.id} $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n 1 | sed -e "s/samtools //g")
        bedtools: \$(bedtools --version | head -n 1 | sed -e "s/bedtools //g")
    END_VERSIONS

    """
}
