process ANNOTATE_CONSENSUS_PEAKS {
    module = [ 'Homer/4.11-Perl-5.30.0' ]

    label "process_high"

    tag "Annotate consensus peaks on $conp"

    input:
    val ( genome )
    path ( gtf )
    path (conp)

    output:
    tuple path(conp), path ( "*.annotation.txt" ), emit: txt
    path ( gtf ), emit: gtf, optional: true
    path ( "*" )

    script:
    def gtf = gtf.baseName != 'dummy_file.txt' ? "$gtf" : ''
    """
    annotate_consensus_peaks.sh $genome $conp $gtf

    """

}
