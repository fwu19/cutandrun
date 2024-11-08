process READS_IN_CONSENSUS_PEAKS {
    module = ['Subread/2.0.0-GCC-8.3.0']

    label "process_high"

    tag "Count reads from ${meta.id} mapped within $conp"

    input:
    tuple val(meta), path(bam), path( conp, stageAs: "conp/*" )

    output:
    tuple val(meta.target), path ( "*.txt" ), emit: count
    path ( "*.fragmentCounts.*" )

    script:
    """
    featureCounts.sh ${meta.id} $bam

    """

}
