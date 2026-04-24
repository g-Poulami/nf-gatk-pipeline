/*
    SAMTOOLS

    SAMTOOLS_SORT    : SAM -> sorted BAM (samtools sort accepts SAM directly)
    SAMTOOLS_INDEX   : sorted BAM -> BAI index
    SAMTOOLS_FLAGSTAT: sorted BAM -> alignment statistics (consumed by MultiQC)
*/

process SAMTOOLS_SORT {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/samtools/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(sam)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: sorted_bam

    script:
    def threads = params.sort_threads ?: task.cpus
    """
    samtools sort \\
        -@ ${threads} \\
        -o ${meta.id}.sorted.bam \\
        ${sam}
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    """
}

process SAMTOOLS_INDEX {
    tag        "${meta.id}"
    label      'process_single'
    publishDir "${params.outdir}/samtools/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bai"), emit: bai

    script:
    """
    samtools index ${bam}
    """

    stub:
    """
    touch ${bam}.bai
    """
}

process SAMTOOLS_FLAGSTAT {
    tag        "${meta.id}"
    label      'process_single'
    publishDir "${params.outdir}/samtools/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.flagstat"), emit: flagstat

    script:
    """
    samtools flagstat ${bam} > ${meta.id}.flagstat
    """

    stub:
    """
    echo "100000 + 0 in total (QC-passed reads + QC-failed reads)" > ${meta.id}.flagstat
    """
}
