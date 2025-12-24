process GENERATE_REPORT {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Make plots of read and peak metrics "

    input:
    path( samplesheet, stageAs: "sample_sheet.csv" )
    path( read_metrics )
    path( "fragment_lengths/*" )
    path( "original_peaks/?/*")
    path( "original_peak_widths/?/*")
    path( "reads_in_peak/?/*" )
    path( "replicated_peaks/*" )
    path( "consensus_peaks/*" )
    path( "consensus_beds/*" )
    path( "consensus_annotation/*" )
    path( "differential_peaks/*" )
    path( "report" )

    output:
    path ( "data.rds" )
    path ( "*.{Rmd,html}" )
    path ( 'versions.yml' ), emit: versions

    script:
    def args = task.ext.args ?: ''
    def suffix = task.ext.suffix ?: ""
    """
    prepare_report.r
    mv report.Rmd 00_analysis_report${suffix}.Rmd
    render_report.r 00_analysis_report${suffix}.Rmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """
}
