process GENERATE_REPORT {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Make plots of read and peak metrics "

    input:
    path( samplesheet, stageAs: "sample_sheet.csv" )
    path( read_metrics )
    path( "fragment_lengths/?/*" )
    path( opeaks, stageAs: "original_peaks/?/*")
    path( "original_peak_widths/?/*")
    path( "reads_in_peak/?/*" )
    path( "replicated_peaks/?/*" )
    path( "consensus_peaks/?/*" )
    path( "consensus_beds/?/*" )
    path( "consensus_annotation/?/*" )
    path( "differential_peaks/?/*" )
    path( "*" )

    output:
    tuple path( "*.{rds,html,Rmd}" )

    script:
    def args = task.ext.args ?: ''
    def suffix = task.ext.suffix ?: ""
    """
    prepare_report.r
    mv report.Rmd 00_analysis_report${suffix}.Rmd
    render_report.r 00_analysis_report${suffix}.Rmd
    """
}
