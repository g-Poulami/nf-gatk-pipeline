/*
    TRIMMOMATIC
    Input : [ meta, [R1, R2] ], adapter_fasta
    Output: paired trimmed reads, unpaired reads, log
*/

process TRIMMOMATIC {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/trimmomatic/${meta.id}", mode: 'copy',
               saveAs: { filename -> filename.endsWith('.log') ? filename : null }

    input:
    tuple val(meta), path(reads)
    path  adapters

    output:
    tuple val(meta), path("*_{1,2}P.fastq.gz"),  emit: trimmed_reads
    tuple val(meta), path("*_{1,2}U.fastq.gz"),  emit: unpaired_reads
    tuple val(meta), path("*.trimmomatic.log"),   emit: log

    script:
    def prefix = meta.id
    """
    trimmomatic PE \\
        -threads ${task.cpus} \\
        -phred33 \\
        ${reads[0]} ${reads[1]} \\
        ${prefix}_1P.fastq.gz ${prefix}_1U.fastq.gz \\
        ${prefix}_2P.fastq.gz ${prefix}_2U.fastq.gz \\
        ILLUMINACLIP:${adapters}:2:30:10 \\
        LEADING:${params.trim_leading} \\
        TRAILING:${params.trim_trailing} \\
        SLIDINGWINDOW:${params.trim_slidingwindow} \\
        MINLEN:${params.trim_minlen} \\
        2> ${prefix}.trimmomatic.log
    """

    stub:
    """
    touch ${meta.id}_1P.fastq.gz ${meta.id}_2P.fastq.gz
    touch ${meta.id}_1U.fastq.gz ${meta.id}_2U.fastq.gz
    touch ${meta.id}.trimmomatic.log
    """
}
