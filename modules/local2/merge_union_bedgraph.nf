
process MERGE_UNION_BEDGRAPH {
    time = '1d'
    cpus = 8
    memory = '48G'

    module = ['fhR/4.1.2-foss-2021b']

    tag "merge union bedgraph files "

    publishDir "${params.outdir}/merged_bedgraph/", mode: 'copy'

    input:
    tuple val(id), path(ubdg)

    output:
    tuple val(id), path("${id}*.bedGraph"), emit: bedgraph
    path("${id}*.bedGraph")

    script:
    def args = task.ext.args ?: ""
    def suffix = task.ext.suffix ?: ""
    """
    merge_union_bedgraph.r $ubdg ${id}${suffix}.bedGraph

    """
}
