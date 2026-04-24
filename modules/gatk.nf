/*
    GATK 4 modules

    GATK_DICT             : CreateSequenceDictionary — required by all GATK tools
    GATK_FAIDX            : samtools faidx — creates .fai index for reference
    GATK_BASERECALIBRATOR : First pass of BQSR — learns error model from known sites
    GATK_APPLYBQSR        : Second pass — applies corrections to BAM
    GATK_HAPLOTYPECALLER  : Per-sample variant calling in gVCF mode
    GATK_GENOTYPEGVCFS    : Joint genotyping across all per-sample gVCFs (optional)
    GATK_VARIANTFILTRATION: Hard filters on SNPs and indels

    Why BQSR?
    ---------
    Illumina assigns Phred quality scores to every base call, but these scores are
    systematically biased by read position, machine cycle, and neighbouring base
    context. BaseRecalibrator learns these biases by comparing the observed
    mismatch rate at known variant sites against the reported quality scores, then
    produces a recalibration table. ApplyBQSR uses that table to correct every
    base's score so HaplotypeCaller receives accurate error probabilities.

    Why gVCF mode in HaplotypeCaller?
    ----------------------------------
    With -ERC GVCF, HaplotypeCaller emits a record for every genomic position,
    not just variant sites. Non-variant positions are compressed into blocks.
    This allows samples to be added to a cohort at any point and re-genotyped
    jointly with GenotypeGVCFs without re-running HaplotypeCaller on old samples.

    Why hard filters instead of VQSR?
    -----------------------------------
    VQSR learns a recalibration model from a truth/training set but requires
    large cohorts (at least 30 WGS samples for SNPs). For smaller projects the
    GATK-recommended approach is hard filters on annotation values such as QD,
    FS, MQ, MQRankSum, and ReadPosRankSum. SNPs and indels have different error
    profiles so they are filtered separately then merged back.
*/

// ── Reference preprocessing ───────────────────────────────────────────────────────

process GATK_DICT {
    tag        "${fasta.baseName}"
    label      'process_single'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "*.dict", emit: dict

    script:
    """
    gatk CreateSequenceDictionary -R ${fasta}
    """

    stub:
    """
    touch ${fasta.baseName}.dict
    """
}

process GATK_FAIDX {
    tag        "${fasta.baseName}"
    label      'process_single'
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    path fasta

    output:
    path "*.fai", emit: fai

    script:
    """
    samtools faidx ${fasta}
    """

    stub:
    """
    touch ${fasta}.fai
    """
}

// ── BQSR ──────────────────────────────────────────────────────────────────────────

process GATK_BASERECALIBRATOR {
    tag        "${meta.id}"
    label      'process_gatk'
    publishDir "${params.outdir}/bqsr/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(index), path(dict), path(fai)
    path  known_sites

    output:
    tuple val(meta), path("*.recal.table"), emit: table

    script:
    def prefix     = meta.id
    def known_args = known_sites instanceof List
        ? known_sites.collect { "--known-sites ${it}" }.join(" \\\n        ")
        : "--known-sites ${known_sites}"
    """
    gatk BaseRecalibrator \\
        -I ${bam} \\
        -R ${fasta} \\
        ${known_args} \\
        -O ${prefix}.recal.table
    """

    stub:
    """
    touch ${meta.id}.recal.table
    """
}

process GATK_APPLYBQSR {
    tag        "${meta.id}"
    label      'process_gatk'
    publishDir "${params.outdir}/bqsr/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(table), path(fasta), path(index), path(dict), path(fai)

    output:
    tuple val(meta), path("*.bqsr.bam"),     emit: bam
    tuple val(meta), path("*.bqsr.bam.bai"), emit: bai

    script:
    def prefix = meta.id
    """
    gatk ApplyBQSR \\
        -I ${bam} \\
        -R ${fasta} \\
        --bqsr-recal-file ${table} \\
        -O ${prefix}.bqsr.bam

    samtools index ${prefix}.bqsr.bam
    """

    stub:
    """
    touch ${meta.id}.bqsr.bam ${meta.id}.bqsr.bam.bai
    """
}

// ── HaplotypeCaller ───────────────────────────────────────────────────────────────

process GATK_HAPLOTYPECALLER {
    tag        "${meta.id}"
    label      'process_gatk'
    publishDir "${params.outdir}/haplotypecaller/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(index), path(dict), path(fai)

    output:
    tuple val(meta), path("*.g.vcf.gz"),     emit: gvcf
    tuple val(meta), path("*.g.vcf.gz.tbi"), emit: tbi

    script:
    def prefix    = meta.id
    def extra     = params.hc_extra_args ?: ''
    def dbsnp_arg = params.dbsnp ? "--dbsnp ${params.dbsnp}" : ''
    """
    gatk HaplotypeCaller \\
        -I ${bam} \\
        -R ${fasta} \\
        -ERC GVCF \\
        ${dbsnp_arg} \\
        ${extra} \\
        -O ${prefix}.g.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.g.vcf.gz ${meta.id}.g.vcf.gz.tbi
    """
}

// ── Joint genotyping (optional) ───────────────────────────────────────────────────

process GATK_GENOTYPEGVCFS {
    label      'process_gatk'
    publishDir "${params.outdir}/genotypegvcfs", mode: 'copy'

    input:
    // collected list of all per-sample gVCFs, then ref tuple
    tuple path(gvcfs), path(fasta), path(index), path(dict), path(fai)

    output:
    path "cohort.vcf.gz",     emit: vcf
    path "cohort.vcf.gz.tbi", emit: tbi

    script:
    def vcf_args = (gvcfs instanceof List ? gvcfs : [gvcfs])
        .findAll { it.name.endsWith('.vcf.gz') }
        .collect { "-V ${it}" }
        .join(" \\\n        ")
    """
    gatk CombineGVCFs \\
        -R ${fasta} \\
        ${vcf_args} \\
        -O combined.g.vcf.gz

    gatk GenotypeGVCFs \\
        -R ${fasta} \\
        -V combined.g.vcf.gz \\
        -O cohort.vcf.gz
    """

    stub:
    """
    touch cohort.vcf.gz cohort.vcf.gz.tbi
    """
}

// ── Variant filtration ────────────────────────────────────────────────────────────

process GATK_VARIANTFILTRATION {
    tag        "${vcf}"
    label      'process_gatk'
    publishDir "${params.outdir}/variantfiltration", mode: 'copy'

    input:
    tuple path(vcf), path(tbi), path(fasta), path(index), path(dict), path(fai)

    output:
    path "*.filtered.vcf.gz",     emit: vcf
    path "*.filtered.vcf.gz.tbi", emit: tbi

    script:
    def prefix    = vcf.simpleName
    def snp_f     = params.snp_filter
    def indel_f   = params.indel_filter
    """
    gatk SelectVariants \\
        -R ${fasta} -V ${vcf} \\
        --select-type-to-include SNP \\
        -O snps.vcf.gz

    gatk VariantFiltration \\
        -R ${fasta} -V snps.vcf.gz \\
        --filter-expression "${snp_f}" \\
        --filter-name "snp_hard_filter" \\
        -O snps.filtered.vcf.gz

    gatk SelectVariants \\
        -R ${fasta} -V ${vcf} \\
        --select-type-to-include INDEL \\
        -O indels.vcf.gz

    gatk VariantFiltration \\
        -R ${fasta} -V indels.vcf.gz \\
        --filter-expression "${indel_f}" \\
        --filter-name "indel_hard_filter" \\
        -O indels.filtered.vcf.gz

    gatk MergeVcfs \\
        -I snps.filtered.vcf.gz \\
        -I indels.filtered.vcf.gz \\
        -O ${prefix}.filtered.vcf.gz
    """

    stub:
    """
    touch ${vcf.simpleName}.filtered.vcf.gz ${vcf.simpleName}.filtered.vcf.gz.tbi
    """
}
