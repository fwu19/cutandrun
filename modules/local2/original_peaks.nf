process ORIGINAL_PEAKS {
    label "process_single"

    tag "Collect metrics of original peaks from ${meta.id}"

    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    input:
    tuple val(meta), path( "peaks/*" )

    output:
    tuple val(meta), path( "*.csv" ), emit: csv
    path ('versions.yml'), emit: versions

    script:
    """
    original_peaks.r ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """

}
