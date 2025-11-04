process GENERATE_REPORT_PROCESS_CONTROLS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Make plots of read and peak metrics "

    input:
    path( samplesheet, stageAs: "sample_sheet.csv" )
    path( read_metrics )
    path( "fragment_lengths/*" )
    path( "original_peaks/*")
    path( "original_peak_widths/*")
    path( "reads_in_peak/*" )
    path( "replicated_peaks/*" )
    path( "*" )
    path( "report/*" )

    output:
    path( "*.{rds,html,Rmd}" ), optional: true
    path ( 'versions.yml' ), emit: versions

    script:
    """
    prepare_report_process_controls.r
    mv report/*.Rmd .
    render_report.r *.Rmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """
}
