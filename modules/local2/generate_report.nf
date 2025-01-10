process GENERATE_REPORT {
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
    path( "consensus_peaks/*" )
    path( "consensus_beds/*" )
    path( "consensus_annotation/*" )
    path( "differential_peaks/*" )
    path( "*" )

    output:
    tuple path( "*.{rds,html,Rmd}" )

    script:
    """
    prepare_report.r
    mv report.Rmd 00_analysis_report.Rmd
    render_report.r 00_analysis_report.Rmd
    """
}
