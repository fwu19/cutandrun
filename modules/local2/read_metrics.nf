process READ_METRICS {

    label "process_single"
    tag "Collect reads QC"

    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    input:
    path ("multiqc_data/*")

    output:
    path( "*.csv" ), emit: csv
    // read_metrics.csv

    script:
    """
    read_metrics.r

    """
    // input: multiqc_data/
}
