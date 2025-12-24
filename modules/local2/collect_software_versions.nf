process COLLECT_SOFTWARE_VERSIONS {
    label 'process_single'

    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    tag "MultiQC on software versions"

    input:
    path (multiqc_config)
    path ('software_versions/*')

    output:
    path ('*multiqc_data/*'), emit: data

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
