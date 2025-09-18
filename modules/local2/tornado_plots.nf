process TORNADO_PLOTS {
    tag "tornado plots"
    label 'process_medium'

    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    input:
    tuple val(target), path("bigwig/*"), path("bed/*")
    path("gene/*")

    output:
    path( "*.pdf" ), optional: true
    path( "*.bed" ), optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    tornado_plots.r target=$target $args
    cp gene/*.bed .
    """
}
