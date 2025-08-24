# ORF15-pipeline

NGS data analysis focused on RPGR-ORF15

## Introduction

Variants in **RPGR** are the leading cause of X-linked retinitis pigmentosa (XLRP). Most pathogenic variants occur in the **ORF15 exon**, a purine-rich, low-complexity region prone to sequencing and alignment errors, complicating accurate variant interpretation. Standard short-read methods such as **WES** and **WGS** often fail to provide sufficient coverage, resulting in unreliable calls or false negatives.

To address this challenge, we developed a refined workflow combining **long-range PCR (LR-PCR)** with optimized bioinformatics. The pipeline applies stringent quality control, alignment filtering, and consensus variant calling, ensuring higher specificity and improved detection of **RPGR-ORF15 variants** for both research and clinical diagnostics.

## Workflow Overview

### 1. Tools & Versions

Tools and reference used:

- fastp - v0.23.2 (quality control)  
- seqtk - (downsampling)  
- BWA-MEM - v0.7.17 (read alignment to hg19)  
- Samclip - (remove reads with clipped bases ≥ 5)  
- SAMtools - v1.13 (filter reads with MAPQ < 60)  
- Picard - v2.25.7 (MarkDuplicates)  
- GATK - v4.2.0.0 (HaplotypeCaller, Mutect2)  
- Python (custom script) - v3.10  
- Reference genome - hg19

### 2. File Structure

```plaintext
your_ORF15_path/

  ORF15_merged_fastq.sh
  sh_all.py
  variant_voting.py
  vcf_1.py
  VCF_processing.py
  config.json
  list.txt              # list of ${sample_name}

  fastq/                # input FASTQ files
  run/                  # auto-generated folder for all analysis outputs
```

### 3. Configuration & Input

- **config.json** → set tool paths here  
- **list.txt** → list of sample names (one per line)  
- **fastq/** → ensure all FASTQ files follow naming format: **${sample_name}-PCR_R1.fastq.gz** & **${sample_name}-PCR_R2.fastq.gz**

### 4. Usage

```bash
# Step 1: Process FASTQ reads
sh ORF15_merged_fastq.sh 

# Step 2: Main analysis pipeline
python3 sh_all.py      # Scripts will automatically be submitted (sbatch) to the HPC cluster for execution

# Step 3: Process VCF file
python3 vcf_1.py 

# Step 4: Generate variant results
python3 VCF_processing.py 

# Step 5: Final consensus calling
python3 variant_voting.py
```

Author: Chien-Yu, Lin
