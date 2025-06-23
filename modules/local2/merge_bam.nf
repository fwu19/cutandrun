
process MERGE_BAM {
    time = '1d'
    cpus = 8
    memory = '48G'
    module = ['SAMtools/1.17-GCC-12.2.0']


    tag "merge bam files for ${id}"

    publishDir "${params.outdir}/merged_bam/", mode: 'copy'

    input:
    tuple val(id), path("input/*")

    output:
    tuple val(id), path("${id}*.bam"), emit: bam
    tuple val(id), path("${id}*.bai"), emit: bai
    path("${id}*.{bam,bai}")

    script:
    def suffix = task.ext.suffix ?: ""
    def args = task.ext.args ?: ""
    """
    samtools merge ${args} -@ ${task.cpus} -f -o merged.bam input/*.bam
    samtools sort -@ ${task.cpus} -o ${id}${suffix}.bam merged.bam
    samtools index ${id}${suffix}.bam

    """
}
