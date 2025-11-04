process DIFFERENTIAL_PEAKS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Call differential peaks on $target"

    input:
    tuple val(target), path ( "counts/*" ), path ( "conp/*" ), path (samplesheet), path (comparison)

    output:
    path ( "*.rds" ), emit: data, optional:true
    path ('versions.yml'), emit: versions
    path ( "*" ), optional:true

    script:
    def args = task.ext.args ?: ''
    """
    differential_peaks.r ss_csv=$samplesheet cmp_file=$comparison tgt=$target $args
    rm -r $samplesheet $comparison counts/ conp/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """

}
