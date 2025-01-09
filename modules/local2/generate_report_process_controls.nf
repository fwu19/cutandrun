process GENERATE_REPORT_PROCESS_CONTROLS {
    module = ['fhR/4.1.2-foss-2021b']

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
    path( rmd )

    output:
    path( "*.{rds,html,Rmd}" ), optional: true

    script:
    """
    prepare_report_process_controls.r
    render_report.r ${rmd}

    """
}
