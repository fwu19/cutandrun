
process UNION_BEDGRAPH {
    time = '1d'
    cpus = 8
    memory = '48G'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'biocontainers/bedtools:2.30.0--hc088bd4_0' }"


    tag "union bedgraph files "

    publishDir "${params.outdir}/merged_bedgraph/", mode: 'copy'

    input:
    tuple val(id), path("input/*")

    output:
    tuple val(id), path("${id}*.unionbdg"), emit: bedgraph

    script:
    def args = task.ext.args ?: ""
    def suffix = task.ext.suffix ?: ""
    """
    bedtools unionbedg $args -i input/* >${id}${suffix}.unionbdg

    """
}
