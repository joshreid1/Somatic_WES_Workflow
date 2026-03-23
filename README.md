# Somatic WES Variant Calling & Annotation Pipeline

A Nextflow (DSL2) pipeline for post-processing, merging, and annotating somatic variants
from [nf-core/sarek](https://nf-co.re/sarek) output. Supports both **tumor-only** and
**paired tumor/normal** WES sample types, using three somatic variant callers (Mutect2,
Strelka2, FreeBayes) with multi-tool concordance scoring.

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Summary](#pipeline-summary)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input](#input)
  - [Parameters](#parameters)
  - [Sample Sheet](#sample-sheet)
  - [Expected nf-core/sarek Directory Structure](#expected-nf-coresarek-directory-structure)
- [Usage](#usage)
- [Pipeline Processes](#pipeline-processes)
- [Output](#output)
- [Filtering Criteria](#filtering-criteria)
- [Notes on Tumor-Only vs Paired Samples](#notes-on-tumor-only-vs-paired-samples)

---

## Overview

This pipeline takes variant call output from [nf-core/sarek](https://nf-co.re/sarek) and:

1. Merges somatic variant calls from **Mutect2**, **Strelka2**, and **FreeBayes**
2. Re-genotypes variants with `bcftools mpileup` for depth/strand validation
3. Annotates variants with **VEP**, **gnomAD v4**, **ClinVar**, **CADD**, and **SpliceAI**
4. Scores multi-tool concordance (`Tool_Count`)
5. Filters variants against a curated somatic gene list and quality thresholds
6. Generates a per-sample Excel summary report and a subset CRAM (BAMlet)

---

## Pipeline Summary

```
nf-core/sarek output
        │
        ▼
 [ProcessSample]          Merge Mutect2 + Strelka2 + FreeBayes VCFs
        │                 Normalise multiallelic sites (bcftools norm)
        │                 Somatic filtering (AD thresholds, sample reheadering)
        ▼
  [Split_VCF]             Split merged VCF for parallelised annotation
        │
        ▼
   [Annotate]             bcftools mpileup re-genotyping
        │                 VEP (v112, GRCh38, cache)
        │                 gnomAD v4.0.0 annotation (vcfanno)
        ▼
 [CADD_Run_Container]     CADD scoring (Singularity container)
        ▼
  [Process_CADD]          Add CADD scores to VCF (vcfanno)
        ▼
    [ClinVar]             ClinVar annotation (vcfanno)
        ▼
  [SpliceAI_Run]          SpliceAI splice-site impact scoring
        ▼
   [Join_VCF]             Rejoin split VCF chunks (snpSift + bcftools sort)
        ▼
  [Check_Tools]           Add per-variant Tool_Count field (vcfanno)
        ▼
 [Filter_Variants]        Gene list + quality filter
        ▼
 [Create_Report]          Excel summary (R script)
        ▼
 [Create_BAMlet]          Subset CRAM over gene regions (samtools)
```

---

## Requirements

### Software / Modules

| Tool | Version | Used in |
|---|---|---|
| Nextflow | ≥ 22.10 | Pipeline orchestration |
| bcftools | 1.20 | VCF processing, filtering, normalisation |
| htslib (bgzip/tabix) | 1.20 | VCF indexing |
| ensembl-vep | 112 | Variant annotation |
| vcfanno | 0.3.2 | gnomAD, ClinVar, CADD annotation |
| snpSift | any | VCF splitting and joining |
| samtools | any | BAMlet creation |
| CADD | via Singularity | CADD scoring |
| SpliceAI | via Singularity | Splice-site scoring |
| R | 4.5.1 | Excel report generation |
| Python | 3.x | Strelka2 VAF calculation script |

### Singularity Containers

- **CADD**: `/vast/projects/reidj-project/containers/cadd-scoring_latest.sif`
- **SpliceAI**: `/stornext/Bioinf/data/lab_bahlo/users/reid.j/analysis/Apptainer/spliceai.img`

### Reference Files

| File | Description |
|---|---|
| `Homo_sapiens_assembly38.fasta` | GRCh38 reference FASTA (with `.fai`) |
| VEP cache (v104, GRCh38) | `/stornext/Bioinf/data/lab_bahlo/ref_db/vep-cache/` |
| gnomAD v4.0.0 TOML | vcfanno config for gnomAD annotation |
| ClinVar TOML | vcfanno config (`clinvar_20250330.toml`) |
| CADD annotations | Mounted into container at `/CADD-scripts/data/annotations` |
| SpliceAI gene annotation | `gencode.v38.annotation.txt` |
| Somatic gene list BED | Custom BED file for filtering (`Somatic_Gene_List_v2025-06.bed`) |

---

## Installation

```bash
git clone <repository-url>
cd <repository-name>
```

Ensure Nextflow is available and all required modules/containers are accessible on your HPC.

---

## Input

### Parameters

| Parameter | Default | Description |
|---|---|---|
| `--output_prefix` | `./somatic_nextflow` | Output directory prefix |
| `--ref_fasta` | *(see script)* | Path to GRCh38 reference FASTA |
| `--sample_info` | *(see script)* | Path to TSV sample sheet |
| `--gene_bed` | *(see script)* | BED file of somatic genes of interest |
| `--spliceai_distance` | `500` | Max distance (bp) from splice site for SpliceAI |

Parameters can be overridden at the command line:

```bash
nextflow run main.nf \
  --sample_info /path/to/samples.tsv \
  --output_prefix ./my_run \
  --gene_bed /path/to/gene_list.bed \
  --spliceai_distance 500
```

### Sample Sheet

The `--sample_info` TSV must have a **header row** with the following columns (tab-separated):

| Column | Description |
|---|---|
| `ID` | Unique sample/pair identifier |
| `BLOOD` | Normal (blood) sample name; use `NA` for tumor-only |
| `BRAIN` | Tumor (brain/tissue) sample name |
| `BAM_PATH` | Root directory containing recalibrated BAM/CRAM files |
| `VCF_PATH` | Root directory of nf-core/sarek variant call output |

**Example — paired tumor/normal:**

```tsv
ID	BLOOD	BRAIN	BAM_PATH	VCF_PATH
SAMPLE001	SAMPLE001_blood	SAMPLE001_brain	/path/to/bam_root	/path/to/sarek_results/variant_calling
```

**Example — tumor-only:**

```tsv
ID	BLOOD	BRAIN	BAM_PATH	VCF_PATH
SAMPLE002	NA	SAMPLE002_brain	/path/to/bam_root	/path/to/sarek_results/variant_calling
```

### Expected nf-core/sarek Directory Structure

The pipeline resolves VCF paths from `VCF_PATH` using the following nf-core/sarek conventions:

**Paired tumor/normal:**
```
{VCF_PATH}/
  mutect2/{ID}/{ID}.mutect2.filtered.vcf.gz
  strelka/{BRAIN}_vs_{BLOOD}/{BRAIN}_vs_{BLOOD}.strelka.somatic_snvs.vcf.gz
  strelka/{BRAIN}_vs_{BLOOD}/{BRAIN}_vs_{BLOOD}.strelka.somatic_indels.vcf.gz
  freebayes/{BRAIN}_vs_{BLOOD}/{BRAIN}_vs_{BLOOD}.freebayes.vcf.gz
```

**Tumor-only:**
```
{VCF_PATH}/
  mutect2/{BRAIN}/{BRAIN}.mutect2.filtered.vcf.gz
  freebayes/{BRAIN}/{BRAIN}.freebayes.vcf.gz
```

BAM/CRAM files are resolved from `{BAM_PATH}/{BRAIN}/{BRAIN}.recal.bam` (or `.recal.cram`).

---

## Usage

```bash
nextflow run main.nf \
  --sample_info samples.tsv \
  --output_prefix ./results_run1 \
  -profile singularity \
  -resume
```

Use `-resume` to restart from the last successful checkpoint after failures.

---

## Pipeline Processes

### `ProcessSample`
- Resolves BAM/CRAM and VCF paths for each sample
- For **paired** samples: adds VAF fields to Strelka2 SNV/indel VCFs using a custom Python script, concatenates SNV+indel VCFs, then reheaders and filters variants requiring `AD_tumor > 1` and `AD_normal < 1`
- For **tumor-only** samples: filters Mutect2 and FreeBayes for `AD[alt] ≥ 1`; Strelka2 is excluded (empty placeholder created)
- Normalises multiallelic variants (`bcftools norm`)
- Tags each variant with its source caller (`Mutect`, `Strelka`, `Freebayes`) in the INFO field
- Outputs a merged VCF combining all callers

### `Split_VCF`
- Splits the merged VCF into chunks using `snpSift split` for parallelised downstream annotation
- Chunk size is adaptive: 500 variants/chunk (>1000 variants), 10/chunk (>10), or 2/chunk otherwise

### `Annotate`
- Re-genotypes all variant positions via `bcftools mpileup` on the recalibrated BAM/CRAM
- Filters for strand-supported alt reads: `ADF[alt] ≥ 1`, `ADR[alt] ≥ 1`, `AD[alt] ≥ 3`
- Runs **Ensembl VEP** (v112, GRCh38): adds consequence, HGVS, SIFT, PolyPhen, AF, AF_1KG, max_AF, variant class, and more
- Annotates with **gnomAD v4.0.0** allele counts (vcfanno, two-pass)

### `CADD_Run_Container`
- Runs **CADD scoring** inside a Singularity container
- Extracts variant positions, strips `chr` prefix, and scores via the CADD Snakemake pipeline

### `Process_CADD`
- Indexes CADD output TSV and adds `CADD_Score` annotations via vcfanno

### `ClinVar`
- Annotates with the latest **ClinVar** VCF via vcfanno

### `SpliceAI_Run`
- Predicts splice-site impact using **SpliceAI** (delta score, distance = `params.spliceai_distance`)
- Runs inside a Singularity container with the GENCODE v38 gene annotation

### `Join_VCF`
- Recombines annotated VCF chunks using `snpSift split -j` and `bcftools sort`

### `Check_Tools`
- Uses vcfanno to annotate each variant with per-caller flags (`Mutect`, `Strelka`, `Freebayes`)
- Computes `Tool_Count` (sum of supporting callers; max = 3 for paired, max = 2 for tumor-only)
- Strelka is excluded from the TOML config for tumor-only samples

### `Filter_Variants`
- Restricts to variants overlapping the somatic gene BED file (`--gene_bed`)
- Applies the following quality filters:

| Filter | Threshold |
|---|---|
| gnomAD v4 allele count | `AC_gnomad_total_4.0 ≤ 1000` |
| Tumor alt read depth | `AD[tumor:alt] ≥ 4` |
| Tumor total read depth | `AD[tumor:ref] + AD[tumor:alt] ≥ 10` |
| Strand balance | `ADF[tumor:alt] ≥ 2` AND `ADR[tumor:alt] ≥ 2` |
| Variant distance bias | `VDB ≥ 0.05` |

### `Create_Report`
- Runs an R script (`filter_somatic_variants.R`) to generate a per-sample **Excel (.xlsx)** summary of filtered variants

### `Create_BAMlet`
- Uses `samtools view` to extract reads overlapping the gene BED file
- Outputs a sorted, indexed **subset CRAM** for IGV inspection

---

## Output

All final outputs are published to `results/`:

```
results/
├── {ID}_somatic_variants_merged.vcf.gz                  # Annotated, caller-tagged VCF
├── {ID}_somatic_variants_merged.vcf.gz.tbi
├── {ID}_somatic_variants_filtered_gene_list.vcf.gz      # Quality-filtered VCF
├── {ID}_somatic_variants_filtered_gene_list.vcf.gz.tbi
├── {ID}_somatic_variants_filtered_summary.xlsx          # Excel report
├── {BRAIN}.subset.cram                                  # BAMlet (gene regions)
└── {BRAIN}.subset.cram.crai
```

---

## Filtering Criteria

The `Filter_Variants` process applies sequential filters. The commented-out lines in the
process (MODERATE/HIGH impact, CADD score) can be re-enabled by uncommenting the relevant
`bcftools view` or `grep` steps in the script block.

---

## Notes on Tumor-Only vs Paired Samples

| Behaviour | Tumor-Only (`BLOOD=NA`) | Paired Tumor/Normal |
|---|---|---|
| Strelka2 | Excluded (empty placeholder) | SNV + indel VCFs merged |
| Somatic filter | `AD[alt] ≥ 1` in tumor | `AD_tumor > 1` AND `AD_normal < 1` |
| `Tool_Count` max | 2 (Mutect2 + FreeBayes) | 3 (Mutect2 + Strelka2 + FreeBayes) |
| Sample reheadering | Not required | Samples renamed to `{ID}_{sample}` |

---

## Authors

Josh Reid — Bahlo Lab, Walter and Eliza Hall Institute of Medical Research
