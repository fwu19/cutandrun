process REPLICATED_PEAKS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Generate replicated peaks from $group"

    input:
    tuple val( group ), val( target ), path ( "peaks/*" )
    val ( min_reps )


    output:
    tuple val(target), path ( "*.csv" ), emit: csv, optional: true
    tuple val(target), path ( "*.rds" ), emit: rds, optional: true
    tuple val(target), path ( "{multiple_replicates,single_replicate}/*.bed" ), emit: bed, optional: true
    path ( "versions.yml" ), emit: versions

    script:
    """
    replicated_peaks.r $group $min_reps

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """

}
