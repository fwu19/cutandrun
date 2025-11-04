process UPDATE_PEAKS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Update peaks"

    input:
    path ( "peaks/*" )
    path ( "peak_annotations/*" )
    path ( "differential_peaks/*" )

    output:
    path ( "*.bed" )
    path ('versions.yml'), emit: versions

    script:
    """
    update_peaks.r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """

}
