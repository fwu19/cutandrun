process AWK {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path( input )

    output:
    tuple val(meta), path("${prefix}.${suffix}"), emit: file
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix   = task.ext.suffix ? "${task.ext.suffix}" : "txt"
    def command  = task.ext.command ?: ''
    def command2 = task.ext.command2 ?: ''

    shell:
    '''
    #!/usr/bin/env bash

    awk -F "\t" '$7 > 2 && $9 > 2' $input $command2 > ${prefix}.${suffix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>/dev/null | head -n 1 | awk '{split(\$0,a,","); print a[1];}' | egrep -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    '''
}
