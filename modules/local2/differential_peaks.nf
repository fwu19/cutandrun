process DIFFERENTIAL_PEAKS {
    //module = ['fhR/4.1.2-foss-2021b']
    container "docker://fwu19/r-libs:4.1.2"

    label "process_single"

    tag "Call differential peaks on $target"

    input:
    tuple val(target), path ("counts/*"), path (samplesheet), path (comparison)

    output:
    path ( "*.dp.rds" ), emit: data, optional:true
    path ('versions.yml'), emit: versions
    path ( "*_peaks" ), optional:true

    script:
    def args = task.ext.args ?: ''
    """
    differential_peaks.r ss_csv=$samplesheet cmp_file=$comparison tgt=$target fdr=${params.fdr} fc=${params.fc} fdr2=${params.fdr2} fc2=${params.fc2} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n 1)
    END_VERSIONS

    """

}
