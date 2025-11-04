process READ_METRICS {

    label "process_single"
    tag "Collect reads QC"

    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    input:
    path ("multiqc_data/*")

    output:
    path( "*.csv" ), emit: csv
    path ('versions.yml'), emit: versions

    script:
    """
    read_metrics.r

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """
    
}
