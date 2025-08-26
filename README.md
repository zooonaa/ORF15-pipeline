# ORF15-pipeline

NGS data analysis focused on RPGR-ORF15

## Introduction

Variants in **RPGR** are the leading cause of X-linked retinitis pigmentosa (XLRP). Most pathogenic variants occur in the **ORF15 region**, a purine-rich, low-complexity region prone to sequencing and alignment errors, complicating accurate variant interpretation. Standard short-read methods such as **WES** and **WGS** often fail to provide sufficient coverage, resulting in unreliable calls or false negatives.

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

#### Scripts & files

- ORF15_merged_fastq.sh: A bash script processing all fastq.gz files. Included QC, merged paired-end data, down-sampling, generated 100 batches fastq file.
- sh_all.py: A python script generated shell scripts and will directly submitted(sbatch) the job to HPC cluster. 
- vcf_1.py: A python script 
- VCF_processing.py
- variant_voting.py
  
- config.json: Please set your tool paths here
- list.txt: List of sample names to be analyzed (each per line) 

### 2. Overall File Structure

```plaintext
your_ORF15_path/
  ORF15_merged_fastq.sh
  sh_all.py
  variant_voting.py
  vcf_1.py
  VCF_processing.py
  config.json
  list.txt              # list of ${sample_name}

  fastq/                # input FASTQ files here: *fastq.gz file under fastq/ is required. ensure all FASTQ files follow naming format: **${sample_name}-PCR_R1.fastq.gz** & **${sample_name}-PCR_R2.fastq.gz**
  run/                  # auto-generated folder for all analysis outputs
```

### 3. Usage

1. Put all scripts under your_ORF15_path/
2. Ensure all FASTQ files follow naming format: **${sample_name}-PCR_R1.fastq.gz** & **${sample_name}-PCR_R2.fastq.gz** and put in the your_ORF15_path/fastq/ folder
3. Make sure config.json and list.txt is filled with correct information

4.
```bash
# Step 1: Process FASTQ reads
sh ORF15_merged_fastq.sh 

# Step 2: Main analysis pipeline
python3 sh_all.py

# Step 3: Process VCF file
python3 vcf_1.py 

# Step 4: Generate variant results
python3 VCF_processing.py 

# Step 5: Final consensus calling
python3 variant_voting.py
```

#### Example of one main script sh_all.py generated:

Sample name = HRD166
Batch = 1 (1 in 100 batch)

```bash
#!/usr/bin/sh
###! You may add the inquiry? the HPC cluster needed in sh_all.py script.

/path/to/BWA/BWA_v0.7.17/bwa mem -t 16 -T 100 -R '@RG\tID:HRD166_PCR_bwamem\tLB:HRD166_PCR_bwamem\tSM:HRD166_PCR_bwamem\tPL:ILLUMINA' /path/to/human_g1k_v37_decoy.fasta /your_ORF15_path/fastq_trimmed/d_single_test100/HRD166-PCR_merge_qc_s1_100000.fastq > /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.sam &&

/path/to/samclip --ref /path/to/human_g1k_v37_decoy.fasta.fai < /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.sam > /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_samclip_.sam &&

java -Xmx80g -jar /path/to/Picard_v2.25.7/picard.jar SortSam INPUT=/your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_samclip_.sam OUTPUT=/your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

/path/to/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools view -h -q 60 -b /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.bam > /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_fq.bam &&

/path/to/SAMTOOLS/SAMTOOLS_v1.13/bin/samtools index /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_fq.bam &&

java -Xmx80g -jar /path/to/Picard_v2.25.7/picard.jar MarkDuplicates INPUT=/your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_fq.bam OUTPUT=/your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.marked.bam METRICS_FILE=/your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_metrics VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

java -Xmx80g -jar /path/to/Picard_v2.25.7/picard.jar SortSam INPUT=/your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.marked.bam OUTPUT=/your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.marked.indexed.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&

java -Xmx80g -jar /path/to/gatk-package-4.2.0.0-local.jar HaplotypeCaller --minimum-mapping-quality 60 -L X:38145300-38145800 -R /path/to/human_g1k_v37_decoy.fasta -I /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.marked.indexed.bam -O /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.haplotype_region.SnpIndel.vcf.gz -bamout /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.haplotype_region.bamout.bam

java -Xmx80g -jar /path/to/gatk-package-4.2.0.0-local.jar Mutect2 --minimum-mapping-quality 60 -L X:38145300-38145800 -R /path/to/human_g1k_v37_decoy.fasta -I /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.marked.indexed.bam -O /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.Mutect2_region.vcf.gz -bamout /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.Mutect2_region.bamout.bam

rm /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.sam
rm /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_samclip_.sam
rm /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem.bam
rm /your_ORF15_path/run/HRD166/s1/HRD166_PCR_bwamem_fq.bam
```

Author: Chien-Yu, Lin
