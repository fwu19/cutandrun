process FILTER_PEAKS {
    tag "filter peaks from $meta.id"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path( input )

    output:
    tuple val(meta), path("*.filtered.*"), emit: file
    path("versions.yml"), emit: versions
    path("*.filtered.*"), optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix   = task.ext.suffix ? "${task.ext.suffix}" : "txt"

    """
    [[ $prefix =~ 'null' ]] || awk -F "\\t" '\$7 > 2 && \$9 > 2' $input > ${prefix}.filtered.${suffix}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk -Wversion 2>/dev/null | head -n 1 | awk '{split(\$0,a,","); print a[1];}' | egrep -o "([0-9]{1,}\\.)+[0-9]{1,}")
    END_VERSIONS
    """
}
