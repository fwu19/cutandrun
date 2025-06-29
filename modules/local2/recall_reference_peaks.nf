process RECALL_REFERENCE_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Compute precision and recall for ${meta.id}"

    input:
    tuple val( meta ), path ( "peaks/*" ), path( "reference" )


    output:
    tuple val(meta), path ( "*.csv" ), emit: csv

    script:
    """
    recall_reference_peaks.r ${meta.id} ${meta.sample_group} ${meta.target}

    """

}
