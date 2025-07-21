process BAM_COVERAGE {
    tag "$id"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(id), path(input), path(input_index)

    output:
    tuple val(id), path("*.bigWig")   , emit: bigwig, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${id}.bigWig"

    """
    bamCoverage \
    --bam $input \
    $args \
    --numberOfProcessors ${task.cpus} \
    --outFileName ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bamCoverage --version | sed -e "s/bamCoverage //g")
    END_VERSIONS
    """
}
