process FRAGMENT_LENGTH {
    label "process_single"

    tag "Compute fragment lengths on ${meta.sample_id}"

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val (meta), path (bam)

    output:
    tuple path("*.fragment_length.txt"), emit: frag_len

    script:
    """
    fragment_length.sh ${meta.id} $bam
    """
}
