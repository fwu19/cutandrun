process SAMPLESHEET_CHECK {
    module = ['fhR/4.1.2-foss-2021b']

    label 'process_single'

    tag "Generate $samplesheet"

    input:
    path ( samplesheet )
    path ( metadata )

    output:
    path '*.csv'        , emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samplesheet_check.r $samplesheet samplesheet.valid.csv $params.use_control $metadata

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}
