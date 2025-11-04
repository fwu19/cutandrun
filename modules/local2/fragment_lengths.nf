process FRAGMENT_LENGTHS {
    label "process_single"

    tag "Compute fragment lengths on ${meta.id}"

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"

    input:
    tuple val (meta), path (bam)

    output:
    tuple val (meta), path("*.fragment_lengths.txt"), emit: txt
    path ( "versions.yml" ), emit: versions

    script:
    """
    fragment_lengths.sh ${meta.id} $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -e "s/samtools //g")
    END_VERSIONS

    """
}
