process READ_METRICS {

    label "process_single"
    tag "Collect reads QC"

    conda "conda-forge::r-base=4.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse:1.2.1' :
        'biocontainers/mulled-v2-03bfeb32fe80910c231f630d4262b83677c8c0f4:f4bb19b68e66de27e4c64306f951d5ff11919931-0' }"

    input:
    path (input)
    path (multiqc_data, stageAs: "multiqc_data/*")
    path (frag_lens, stageAs: "fragment_length/*")

    output:
    path( "*.{csv,rds}" ), emit: metrics


    script:
    """
    read_metrics.r $input

    """
}
