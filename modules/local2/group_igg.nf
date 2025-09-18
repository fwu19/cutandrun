process GROUP_IGG {

    label 'process_single'

    tag "Group IgG"

    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    input:
    path ( "samplesheet.valid.csv" )
    val (igg_group)

    output:
    path ('samplesheet.valid.combine_igg.csv'), emit: csv
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    group_igg.r samplesheet.valid.csv samplesheet.valid.combine_igg.csv $igg_group

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS
    """
}

