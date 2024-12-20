process REPLICATED_HICONF_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Generate replicated peaks from ${meta.id}"

    input:
    tuple val( meta ), path ( "peaks/*" ), path( "*" )


    output:
    tuple val(meta), path ( "*.csv" ), emit: csv

    script:
    """
    replicated_high_confidence_peaks.r ${meta.id} ${meta.sample_group} ${meta.target}

    """

}
