process DIFFERENTIAL_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Call differential peaks on $target"

    input:
    tuple val(target), path ( "counts/*" ), path ( "conp/*" ), path (samplesheet), path (comparison)

    output:
    path ( "*.rds" ), emit: data, optional:true
    path ( "*" ), optional:true

    script:
    def args = task.ext.args ?: ''
    """
    differential_peaks.r ss_csv=$samplesheet cmp_file=$comparison tgt=$target $args
    rm -r $samplesheet $comparison counts/ conp/
    """

    // input files:
    // rm -r $samplesheet $comparison counts/ conp/
}
