process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_single'

    conda "conda-forge::r-base=4.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' :
        'biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0' }"

    input:
    path samplesheet

    output:
    path '*.csv'        , emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_samplesheet.r $samplesheet samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}
