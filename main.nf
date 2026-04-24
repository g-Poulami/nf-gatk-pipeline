#!/usr/bin/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-gatk-pipeline
    BWA-MEM2 -> SAMtools -> Picard MarkDuplicates -> GATK BQSR -> HaplotypeCaller
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW      } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED  } from './modules/fastqc'
include { TRIMMOMATIC               } from './modules/trimmomatic'
include { BWAMEM2_INDEX             } from './modules/bwamem2'
include { BWAMEM2_MEM               } from './modules/bwamem2'
include { SAMTOOLS_SORT             } from './modules/samtools'
include { SAMTOOLS_INDEX            } from './modules/samtools'
include { SAMTOOLS_FLAGSTAT         } from './modules/samtools'
include { PICARD_MARKDUPLICATES     } from './modules/picard'
include { GATK_DICT                 } from './modules/gatk'
include { GATK_FAIDX                } from './modules/gatk'
include { GATK_BASERECALIBRATOR     } from './modules/gatk'
include { GATK_APPLYBQSR            } from './modules/gatk'
include { GATK_HAPLOTYPECALLER      } from './modules/gatk'
include { GATK_GENOTYPEGVCFS        } from './modules/gatk'
include { GATK_VARIANTFILTRATION    } from './modules/gatk'
include { MULTIQC                   } from './modules/multiqc'

log.info """
    nf-gatk-pipeline  v${workflow.manifest.version}
    ================================================
    reads         : ${params.reads}
    genome        : ${params.genome}
    known_sites   : ${params.known_sites ?: 'none (BQSR disabled)'}
    outdir        : ${params.outdir}
    run_bqsr      : ${params.run_bqsr}
    joint_genotype: ${params.joint_genotype}
    """.stripIndent()

def validateParams() {
    def errors = []
    if (!params.reads)  errors << "  --reads is required"
    if (!params.genome) errors << "  --genome is required"
    if (params.run_bqsr && !params.known_sites) {
        errors << "  --known_sites is required when --run_bqsr is true"
    }
    if (errors) {
        log.error "Missing required parameters:\n" + errors.join("\n")
        System.exit(1)
    }
}

workflow {

    validateParams()

    // -----------------------------------------------------------------------
    // Input channels
    // -----------------------------------------------------------------------

    // [ meta, [R1, R2] ]
    ch_reads = Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .map { id, files ->
            [ [id: id, single_end: false, platform: params.sequencing_platform], files ]
        }

    ch_genome = Channel
        .fromPath(params.genome, checkIfExists: true)
        .first()

    ch_adapters = params.adapters
        ? Channel.fromPath(params.adapters, checkIfExists: true).first()
        : Channel.fromPath("${projectDir}/assets/adapters.fa").first()

    // -----------------------------------------------------------------------
    // QC and trimming
    // -----------------------------------------------------------------------

    FASTQC_RAW(ch_reads)
    TRIMMOMATIC(ch_reads, ch_adapters)
    FASTQC_TRIMMED(TRIMMOMATIC.out.trimmed_reads)

    // -----------------------------------------------------------------------
    // Reference indexing — runs once, shared across all samples
    // -----------------------------------------------------------------------

    BWAMEM2_INDEX(ch_genome)
    GATK_DICT(ch_genome)
    GATK_FAIDX(ch_genome)

    // Pack reference + all index files into a single reusable channel.
    // .first() makes it a value channel so it can be combined repeatedly.
    ch_ref = BWAMEM2_INDEX.out.index
        .combine(GATK_DICT.out.dict)
        .combine(GATK_FAIDX.out.fai)
        .first()

    // -----------------------------------------------------------------------
    // Alignment
    // BWA-MEM2 writes SAM; SAMTOOLS_SORT converts and sorts in its own
    // container so no cross-container pipe is required.
    // -----------------------------------------------------------------------

    ch_align_input = TRIMMOMATIC.out.trimmed_reads
        .combine(ch_ref)

    BWAMEM2_MEM(ch_align_input)
    SAMTOOLS_SORT(BWAMEM2_MEM.out.sam)
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.sorted_bam)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT.out.sorted_bam)

    // [ meta, bam, bai ]
    ch_sorted = SAMTOOLS_SORT.out.sorted_bam
        .join(SAMTOOLS_INDEX.out.bai)

    // -----------------------------------------------------------------------
    // Mark duplicates
    // -----------------------------------------------------------------------

    PICARD_MARKDUPLICATES(ch_sorted)

    // [ meta, bam, bai ] after deduplication
    ch_dedup = PICARD_MARKDUPLICATES.out.bam
        .join(PICARD_MARKDUPLICATES.out.bai)

    // -----------------------------------------------------------------------
    // BQSR (optional — requires known variant sites)
    // -----------------------------------------------------------------------

    if (params.run_bqsr) {
        ch_known = Channel
            .fromPath(params.known_sites, checkIfExists: true)
            .collect()

        ch_bqsr_input = ch_dedup.combine(ch_ref)

        GATK_BASERECALIBRATOR(ch_bqsr_input, ch_known)

        ch_apply_input = ch_dedup
            .join(GATK_BASERECALIBRATOR.out.table)
            .combine(ch_ref)

        GATK_APPLYBQSR(ch_apply_input)

        ch_final = GATK_APPLYBQSR.out.bam
            .join(GATK_APPLYBQSR.out.bai)
    } else {
        ch_final = ch_dedup
    }

    // -----------------------------------------------------------------------
    // Variant calling
    // -----------------------------------------------------------------------

    ch_hc_input = ch_final.combine(ch_ref)

    GATK_HAPLOTYPECALLER(ch_hc_input)

    // -----------------------------------------------------------------------
    // Joint genotyping (optional)
    // -----------------------------------------------------------------------

    if (params.joint_genotype) {
        ch_gvcf_input = GATK_HAPLOTYPECALLER.out.gvcf
            .collect()
            .combine(ch_ref)

        GATK_GENOTYPEGVCFS(ch_gvcf_input)
        ch_vcf = GATK_GENOTYPEGVCFS.out.vcf
    } else {
        ch_vcf = GATK_HAPLOTYPECALLER.out.gvcf
    }

    // -----------------------------------------------------------------------
    // Hard variant filtration
    // -----------------------------------------------------------------------

    ch_filter_input = ch_vcf.combine(ch_ref)

    GATK_VARIANTFILTRATION(ch_filter_input)

    // -----------------------------------------------------------------------
    // MultiQC
    // -----------------------------------------------------------------------

    if (params.run_multiqc) {
        ch_qc = Channel.empty()
            .mix(FASTQC_RAW.out.zip)
            .mix(FASTQC_TRIMMED.out.zip)
            .mix(TRIMMOMATIC.out.log)
            .mix(SAMTOOLS_FLAGSTAT.out.flagstat)
            .mix(PICARD_MARKDUPLICATES.out.metrics)
            .collect()

        MULTIQC(ch_qc)
    }
}

workflow.onComplete {
    def status = workflow.success ? "SUCCESS" : "FAILED"
    log.info """
    Pipeline ${status}
    Completed : ${workflow.complete}
    Duration  : ${workflow.duration}
    Output    : ${params.outdir}
    """.stripIndent()
}
