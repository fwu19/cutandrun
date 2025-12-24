process ANNOTATE_CONSENSUS_PEAKS {
    module = [ 'Homer/4.11-Perl-5.30.0' ]

    label "process_high"

    tag "Annotate consensus peaks on $conp"

    input:
    path ( fasta )
    path ( gtf )
    path (conp)

    output:
    tuple path(conp), path ( "*.annotation.txt" ), emit: txt, optional: true
    path ( gtf ), emit: gtf, optional: true
    path ( "versions.yml" ), emit: versions

    script:
    def gtf = gtf.baseName != 'dummy_file.txt' ? "$gtf" : ''
    """
    annotate_consensus_peaks.sh $fasta $gtf $conp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version | head -n 1 | sed -e "s/.*\\( //g; s/\\).*//g")
    END_VERSIONS

    """

}
