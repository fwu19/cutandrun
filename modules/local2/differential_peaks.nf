process DIFFERENTIAL_PEAKS {
    module = ['fhR/4.1.2-foss-2021b']

    label "process_single"

    tag "Call differential peaks on $target"

    input:
    path (samplesheet)
    path (comparison)
    tuple val(target), path ( "counts/*" ), path ( "conp/*" )

    output:
    path ( "*.{rds,txt}" ), emit: data
    path ( "*" )

    script:
    """
    differential_peaks.r $samplesheet $comparison $target
    rm -r $samplesheet $comparison counts/ conp/
    """

    // input files:
    // rm -r $samplesheet $comparison counts/ conp/
}
