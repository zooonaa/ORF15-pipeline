import os, sys
import json

#=================================basic========================================#
script_dir = os.path.dirname(os.path.abspath(__file__))
work_dir   = os.path.join(script_dir, "run")
os.makedirs(work_dir, exist_ok=True)  
# os.chdir(work_dir)

SampleList = os.path.join(script_dir, "list.txt")

fastq_dir = os.path.join(script_dir, "fastq_trimmed")  
fastq_trim_d = os.path.join(fastq_dir, "d_single_test100")

data = "_PCR"
tool = "_all"
align = "_bwamem"
threading = "16"

#=======================software direction=====================================#
with open(os.path.join(script_dir, "config.json")) as f:
    config = json.load(f)

bwa       = os.path.join(config["bwa"])
picard    = os.path.join(config["picard"])
GATK4     = os.path.join(config["GATK4"])
samtools  = os.path.join(config["samtools"])
ref_ucsc  = os.path.join(config["ref_ucsc"])
ref_fai   = os.path.join(config["ref_fai"])
samclip   = os.path.join(config["samclip"])

jar = 'java -Xmx80g -jar'  # 40g/80g
python3 = "/bin/python3.6"

#=============================================================#

with open(SampleList, 'r') as file:
    S_List = [line.strip() for line in file]

node_range = range(1, 101)

for sample in S_List:
    sample_dir = os.path.join(work_dir, sample)
    os.makedirs(sample_dir, exist_ok=True)

    for node in node_range:
        status = f"s{node}"
        status_dir = os.path.join(sample_dir, status)
        os.makedirs(status_dir, exist_ok=True)

        fastq_r1 = os.path.join(fastq_trim_d, f"{sample}-PCR_merge_qc_{status}_100000.fastq")

        sh_path = os.path.join(sample_dir, f"{sample}_{status}.sh")
        with open(sh_path, 'w', encoding='utf-8') as file1:
            file1.write('#!/usr/bin/sh\n')

            sam_path        = os.path.join(status_dir, f"{sample}{data}{align}.sam")
            samclip_path    = os.path.join(status_dir, f"{sample}{data}{align}_samclip_.sam")
            bam_path        = os.path.join(status_dir, f"{sample}{data}{align}.bam")
            bam_fq_path     = os.path.join(status_dir, f"{sample}{data}{align}_fq.bam")
            bam_marked_path = os.path.join(status_dir, f"{sample}{data}{align}.marked.bam")
            bam_sorted_idx  = os.path.join(status_dir, f"{sample}{data}{align}.marked.indexed.bam")
            hc_vcf          = os.path.join(status_dir, f"{sample}{data}{align}.haplotype_region.SnpIndel.vcf.gz")
            hc_bamout       = os.path.join(status_dir, f"{sample}{data}{align}.haplotype_region.bamout.bam")
            mut_vcf         = os.path.join(status_dir, f"{sample}{data}{align}.Mutect2_region.vcf.gz")
            mut_bamout      = os.path.join(status_dir, f"{sample}{data}{align}.Mutect2_region.bamout.bam")
            metrics_file    = os.path.join(status_dir, f"{sample}{data}{align}_metrics")

            # bwa mem
            file1.write(
                f"{bwa} mem -t {threading} -T 100 -R "
                f"'@RG\\tID:{sample}{data}{align}\\tLB:{sample}{data}{align}\\tSM:{sample}{data}{align}\\tPL:ILLUMINA' "
                f"{ref_ucsc} {fastq_r1} > {sam_path} &&\n\n"
            )

            # samclip
            file1.write(f"{samclip} --ref {ref_fai} < {sam_path} > {samclip_path} &&\n\n")

            # picard SortSam
            file1.write(
                f"{jar} {picard} SortSam INPUT={samclip_path} OUTPUT={bam_path} "
                f"SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&\n\n"
            )

            # samtools filter & index
            file1.write(f"{samtools} view -h -q 60 -b {bam_path} > {bam_fq_path} &&\n\n")
            file1.write(f"{samtools} index {bam_fq_path} &&\n\n")

            # picard MarkDuplicates
            file1.write(
                f"{jar} {picard} MarkDuplicates INPUT={bam_fq_path} OUTPUT={bam_marked_path} "
                f"METRICS_FILE={metrics_file} VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&\n\n"
            )

            # picard SortSam (again)
            file1.write(
                f"{jar} {picard} SortSam INPUT={bam_marked_path} OUTPUT={bam_sorted_idx} "
                f"SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true &&\n\n"
            )

            # GATK HaplotypeCaller
            file1.write(
                f"{jar} {GATK4} HaplotypeCaller --minimum-mapping-quality 60 -L X:38145300-38145800 "
                f"-R {ref_ucsc} -I {bam_sorted_idx} -O {hc_vcf} -bamout {hc_bamout}\n\n"
            )

            # GATK Mutect2
            file1.write(
                f"{jar} {GATK4} Mutect2 --minimum-mapping-quality 60 -L X:38145300-38145800 "
                f"-R {ref_ucsc} -I {bam_sorted_idx} -O {mut_vcf} -bamout {mut_bamout}\n\n"
            )

            file1.write(f"rm {sam_path}\n")
            file1.write(f"rm {samclip_path}\n")
            file1.write(f"rm {bam_path}\n")
            file1.write(f"rm {bam_fq_path}\n")

import time
import glob

sh_files = sorted(glob.glob(os.path.join(work_dir, "*/*.sh")))

for sh in sh_files:
    print(f"Submitting {sh} ...")
    os.system(f"sbatch {sh}")
    time.sleep(0.8) 
