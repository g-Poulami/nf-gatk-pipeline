/*
    FASTQC
    Input : [ meta, [R1, R2] ]
    Output: html reports and zip archives
*/

process FASTQC {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/fastqc/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip"),  emit: zip

    script:
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir  . \\
        ${reads}
    """

    stub:
    """
    touch ${meta.id}_R1_fastqc.html ${meta.id}_R1_fastqc.zip
    touch ${meta.id}_R2_fastqc.html ${meta.id}_R2_fastqc.zip
    """
}
