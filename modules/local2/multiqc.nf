process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    tag "MultiQC on all samples"

    input:
    path multiqc_config
    path ('bowtie2/*')
    path ('bowtie2_spikein/*')
    path ('samtools/stats/*')
    path ('samtools/flagstat/*')
    path ('samtools/idxstats/*')
    path ('picard/markduplicates/*')

    output:
    path ('*multiqc*')
    path ('*multiqc_data/*'), emit: data
    path ('*multiqc_report.html'), emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    multiqc -o ./ -f ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
