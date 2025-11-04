process ORIGINAL_PEAK_WIDTHS {
    label "process_single"

    tag "Collect peak widths from ${meta.id}"

    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    input:
    tuple val(meta), path( "peaks/*" )

    output:
    tuple val(meta), path( "*.csv" ), emit: csv
    path ('versions.yml'), emit: versions

    script:
    """
    original_peak_widths.r ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """

}
