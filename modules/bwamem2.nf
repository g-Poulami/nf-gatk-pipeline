/*
    BWA-MEM2

    BWAMEM2_INDEX
    -------------
    Input : reference FASTA
    Output: tuple(fasta, [index_files])
            The FASTA is passed through so downstream processes receive
            the reference and its index together in a single tuple.

    BWAMEM2_MEM
    -----------
    Input : [ meta, [R1, R2], fasta, [index_files], dict, fai ]
    Output: [ meta, sam ]

    This process writes plain SAM, not BAM. Piping bwa-mem2 directly
    into samtools view would require samtools inside the BWA-MEM2
    container, which the biocontainers image does not provide. Keeping
    the output as SAM means SAMTOOLS_SORT (in its own container) handles
    conversion and sorting in a single step. samtools sort accepts SAM
    input natively.
*/

process BWAMEM2_INDEX {
    tag        "${fasta.baseName}"
    label      'process_high'
    publishDir "${params.outdir}/bwamem2_index", mode: 'copy'

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta}.*"), emit: index

    script:
    """
    bwa-mem2 index ${fasta}
    """

    stub:
    """
    touch ${fasta}.0123 ${fasta}.amb ${fasta}.ann
    touch ${fasta}.bwt.2bit.64 ${fasta}.pac
    """
}

process BWAMEM2_MEM {
    tag        "${meta.id}"
    label      'process_high'
    publishDir "${params.outdir}/bwamem2/${meta.id}", mode: 'copy'

    input:
    // Channel layout from main.nf:
    // TRIMMOMATIC.out.trimmed_reads .combine( ch_ref )
    // ch_ref = tuple(fasta, index_files, dict, fai)
    tuple val(meta), path(reads), path(fasta), path(index), path(dict), path(fai)

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
    def prefix = meta.id
    def extra  = params.bwamem2_extra_args ?: ''
    def rg     = "@RG\\tID:${prefix}\\tSM:${prefix}\\tPL:${meta.platform}\\tLB:${prefix}\\tPU:${prefix}"
    """
    bwa-mem2 mem \\
        -t ${task.cpus} \\
        -R "${rg}" \\
        ${extra} \\
        ${fasta} \\
        ${reads[0]} ${reads[1]} \\
        > ${prefix}.sam
    """

    stub:
    """
    touch ${meta.id}.sam
    """
}
