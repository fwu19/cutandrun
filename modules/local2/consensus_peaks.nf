process CONSENSUS_PEAKS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Generate consensus peaks from $target"

    input:
    tuple val (target), path ( "peaks/*" ), path (samplesheet)

    output:
    tuple val(target), path ( "*.csv" ), emit: csv
    tuple val(target), path ( "*.rds" ), emit: rds
    tuple val(target), path ( "*.bed" ), emit: bed
    path ('versions.yml'), emit: versions

    script:
    """
    consensus_peaks.r $samplesheet $target

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featureCounts: \$( featureCounts -v | sed -e "s/featureCounts //g" )
    END_VERSIONS

    """

}
