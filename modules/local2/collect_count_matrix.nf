process COLLECT_COUNT_MATRIX{
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Collect read count matrix on $target"

    input:
    tuple val(target), path ( "counts/*" ), path ( "conp/*" ), path (samplesheet)

    output:
    path ('versions.yml'), emit: versions
    tuple val(target), path ( "*.raw_counts.txt" ), emit: cts, optional: true

    script:
    def args = task.ext.args ?: ''
    """
    collect_count_matrix.r ss_csv=$samplesheet tgt=$target $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """

}
