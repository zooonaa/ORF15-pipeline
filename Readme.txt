Author: Chien-Yu, Lin

Step1: ORF15_merged_fastq.sh (process fastq reads)
Step2: sh_all.py (Run script)
Step3: vcf_1.py (process vcf file)
Step4: VCF_processing.py (generate variant result)
Step5: variant_voting.py (select final decision for all samples)

config1.json: put your tool path here
-ORF15
	-all scripts
	-fastq
	-run