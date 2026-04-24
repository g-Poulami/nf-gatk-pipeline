# nf-gatk-pipeline

BWA-MEM2 -> SAMtools -> Picard -> GATK BQSR -> HaplotypeCaller variant calling pipeline in Nextflow DSL2.

[![CI](https://github.com/g-Poulami/nf-gatk-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/g-Poulami/nf-gatk-pipeline/actions/workflows/ci.yml)
![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen)
![GATK 4.5](https://img.shields.io/badge/GATK-4.5-orange)
![License](https://img.shields.io/badge/license-MIT-blue)

---

## Pipeline steps

```
paired FASTQ
     |
     v
 FastQC (raw)
     |
     v
 Trimmomatic                adapter removal, quality trimming
     |
     v
 FastQC (trimmed)
     |
     v
 BWA-MEM2 index (once)
     |
 BWA-MEM2 mem              align, embed @RG tag, write SAM
     |
     v
 SAMtools sort              SAM -> sorted BAM
     |
     v
 SAMtools index             sorted BAM -> BAI
     |
 SAMtools flagstat          alignment statistics
     |
     v
 Picard MarkDuplicates      mark PCR duplicates
     |
     v
 GATK BaseRecalibrator  <-- known sites VCF (dbSNP, Mills)
     |
 GATK ApplyBQSR             recalibrated BAM
     |
     v
 GATK HaplotypeCaller       per-sample gVCF (-ERC GVCF)
     |
     +-- (--joint_genotype) --> GATK GenotypeGVCFs --> cohort VCF
     |
     v
 GATK VariantFiltration     hard filters: SNPs and indels separately
     |
     v
 Filtered VCF
     |
     v
 MultiQC                    aggregated HTML report
```

### Outputs

| Directory | Contents |
|-----------|----------|
| `results/fastqc/` | FastQC HTML and ZIP (raw and trimmed) |
| `results/trimmomatic/` | Trimmomatic log files |
| `results/bwamem2/` | SAM files |
| `results/samtools/` | Sorted BAM, BAI index, flagstat |
| `results/picard/` | Deduplicated BAM and duplication metrics |
| `results/bqsr/` | Recalibration tables and recalibrated BAM |
| `results/haplotypecaller/` | Per-sample gVCF and TBI index |
| `results/genotypegvcfs/` | Joint-genotyped cohort VCF (if enabled) |
| `results/variantfiltration/` | Filtered VCF |
| `results/multiqc/` | `multiqc_report.html` |
| `results/pipeline_info/` | Timeline, resource report, execution DAG |

---

## Quick start

### Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

### Stub run — no tools required

```bash
git clone https://github.com/g-Poulami/nf-gatk-pipeline.git
cd nf-gatk-pipeline
python3 test/generate_test_data.py
nextflow run main.nf -profile test -stub-run
```

### Run with Docker

```bash
nextflow run main.nf \
  -profile docker \
  --reads        'data/*_R{1,2}.fastq.gz' \
  --genome       'ref/hg38.fa' \
  --known_sites  'ref/dbsnp_146.hg38.vcf.gz' \
  --outdir       results
```

### Skip BQSR (non-human organisms or no known sites available)

```bash
nextflow run main.nf \
  -profile docker \
  --reads    'data/*_R{1,2}.fastq.gz' \
  --genome   'ref/genome.fa' \
  --run_bqsr false \
  --outdir   results
```

### Joint genotyping across multiple samples

```bash
nextflow run main.nf \
  -profile docker \
  --reads          'data/*_R{1,2}.fastq.gz' \
  --genome         'ref/hg38.fa' \
  --known_sites    'ref/dbsnp_146.hg38.vcf.gz' \
  --joint_genotype true \
  --outdir         results
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--reads` | required | Glob pattern for paired FASTQ files |
| `--genome` | required | Reference genome FASTA |
| `--known_sites` | null | VCF(s) for BQSR (required when `--run_bqsr true`) |
| `--dbsnp` | null | dbSNP VCF for HaplotypeCaller annotation |
| `--adapters` | `assets/adapters.fa` | Adapter sequences for Trimmomatic |
| `--outdir` | `results` | Output directory |
| `--run_bqsr` | `true` | Run BQSR |
| `--remove_duplicates` | `false` | Remove duplicates instead of marking |
| `--joint_genotype` | `false` | Run GenotypeGVCFs across all samples |
| `--snp_filter` | GATK defaults | Hard filter expression for SNPs |
| `--indel_filter` | GATK defaults | Hard filter expression for indels |
| `--run_multiqc` | `true` | Aggregate QC reports |

---

## Reference files for GRCh38

| File | Source |
|------|--------|
| `hg38.fa` | UCSC hgdownload |
| `dbsnp_146.hg38.vcf.gz` | GATK Resource Bundle |
| `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | GATK Resource Bundle |

GATK Resource Bundle: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811

---

## Profiles

| Profile | Description |
|---------|-------------|
| `local` | Run locally without containers |
| `docker` | Pull from quay.io/biocontainers and broadinstitute/gatk |
| `singularity` | Use Singularity images (recommended for HPC) |
| `conda` | Per-process conda environments |
| `slurm` | Submit to a SLURM cluster |
| `test` | Synthetic data for CI (`-stub-run`) |

Combine profiles:

```bash
nextflow run main.nf -profile singularity,slurm --reads ...
```

---

## Design notes

**Why SAM output from BWA-MEM2?**
The biocontainers BWA-MEM2 image does not include SAMtools. Piping `bwa-mem2 mem`
directly into `samtools view` in the same script block would fail at runtime
because samtools is not present in that container. Writing SAM and converting in
`SAMTOOLS_SORT` keeps each process entirely within its own container. SAMtools sort
accepts SAM input natively so no extra conversion step is needed.

**Why gVCF mode?**
Running HaplotypeCaller with `-ERC GVCF` stores per-position evidence for every
sample. New samples can be added to the cohort and re-genotyped with
`GenotypeGVCFs` at any time without re-running HaplotypeCaller on existing samples.
This is the standard approach for any project that may grow.

**Hard filters vs VQSR**
VQSR requires at least 30 WGS samples to build a reliable recalibration model.
For smaller cohorts or non-human organisms the GATK-recommended alternative is
hard filters on annotation values. The default thresholds in this pipeline match
GATK's own published recommendations.

---

## Resuming after failure

```bash
nextflow run main.nf -resume [other params]
```

---

## Project structure

```
nf-gatk-pipeline/
├── main.nf
├── nextflow.config
├── assets/
│   └── adapters.fa
├── modules/
│   ├── fastqc.nf
│   ├── trimmomatic.nf
│   ├── bwamem2.nf
│   ├── samtools.nf
│   ├── picard.nf
│   ├── gatk.nf
│   └── multiqc.nf
├── test/
│   ├── generate_test_data.py
│   └── data/
└── .github/
    └── workflows/
        └── ci.yml
```

---

## License

MIT

## Results

Pipeline run on SRR062634 (100,000 reads, chr22 reference). BQSR skipped for test run.

| Metric | Value |
|--------|-------|
| Total reads | 195,301 |
| Mapped reads | 55,273 |
| Mapping rate | 27.7% |

Reports:
- [Pipeline execution report](https://g-Poulami.github.io/nf-gatk-pipeline/report.html)
- [FastQC R1](https://g-Poulami.github.io/nf-gatk-pipeline/fastqc/test_R1_fastqc.html)
- [FastQC R2](https://g-Poulami.github.io/nf-gatk-pipeline/fastqc/test_R2_fastqc.html)
