process BIGWIG_AVERAGE {
    tag "$group"
    label 'process_medium'

    conda "bioconda::deeptools=3.5.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/deeptools:3.5.6--pyhdfd78af_0':
        'biocontainers/deeptools:3.5.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path("input/*")

    output:
    tuple val(meta), path("*.bigWig")   , emit: bigwig, optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta[0]}_${meta[1]}.bigWig"

    """
    bigwigAverage \
    -b input/* \
    $args \
    --numberOfProcessors ${task.cpus} \
    --outFileName ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deeptools: \$(bigwigAverage --version | sed -e "s/bigwigAverage //g")
    END_VERSIONS
    """
}
