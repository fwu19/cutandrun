process RECALL_REFERENCE_PEAKS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Compute precision and recall for ${meta.id}"

    input:
    tuple val( meta ), path ( "peaks/*" ), path( "reference" )


    output:
    tuple val(meta), path ( "*.csv" ), emit: csv
    path ('versions.yml'), emit: versions

    script:
    """
    recall_reference_peaks.r ${meta.id} ${meta.sample_group} ${meta.target}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featureCounts: \$( featureCounts -v | sed -e "s/featureCounts //g" )
    END_VERSIONS

    """

}
