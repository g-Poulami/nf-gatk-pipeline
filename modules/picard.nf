/*
    PICARD MARKDUPLICATES

    PCR duplicates are reads derived from the same original DNA fragment that
    have been amplified multiple times. If left unmarked they inflate allele
    counts and can produce false-positive variant calls. Picard identifies
    duplicates by comparing the outer coordinates of paired reads. It either
    marks them with the SAM FLAG 0x400 (default) or removes them entirely
    depending on params.remove_duplicates.

    Input : [ meta, bam, bai ]
    Output: marked/deduped BAM, BAI index, duplication metrics file
*/

process PICARD_MARKDUPLICATES {
    tag        "${meta.id}"
    label      'process_gatk'
    publishDir "${params.outdir}/picard/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.markdup.bam"),           emit: bam
    tuple val(meta), path("*.markdup.bam.bai"),       emit: bai
    tuple val(meta), path("*.MarkDuplicates.metrics"), emit: metrics

    script:
    def prefix     = meta.id
    def remove_dup = params.remove_duplicates ? 'true' : 'false'
    """
    picard MarkDuplicates \\
        INPUT=${bam} \\
        OUTPUT=${prefix}.markdup.bam \\
        METRICS_FILE=${prefix}.MarkDuplicates.metrics \\
        REMOVE_DUPLICATES=${remove_dup} \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=.

    mv ${prefix}.markdup.bai ${prefix}.markdup.bam.bai
    """

    stub:
    """
    touch ${meta.id}.markdup.bam ${meta.id}.markdup.bam.bai
    echo "LIBRARY\t0.01" > ${meta.id}.MarkDuplicates.metrics
    """
}
