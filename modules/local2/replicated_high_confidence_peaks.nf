process REPLICATED_HICONF_PEAKS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Generate replicated peaks from ${meta.id}"

    input:
    tuple val( meta ), path ( "peaks/*" ), path( "*" )


    output:
    tuple val(meta), path ( "*.csv" ), emit: csv
    path ('versions.yml'), emit: versions

    script:
    """
    replicated_high_confidence_peaks.r ${meta.id} ${meta.sample_group} ${meta.target}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        featureCounts: \$( featureCounts -v | sed -e "s/featureCounts //g" )
    END_VERSIONS

    """

}
