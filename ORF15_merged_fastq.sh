#!/usr/bin/sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/config.json"
if ! command -v jq &> /dev/null; then
    echo "jq not found. Please install jq to use this script."
    exit 1
fi
fastp=$(jq -r '.fastp' "$CONFIG_FILE")
seqtk=$(jq -r '.seqtk' "$CONFIG_FILE")


SampleList="${SCRIPT_DIR}/list.txt"
fastq_path="${SCRIPT_DIR}/fastq"

fastq_raw_a=${fastq_path}/a_rawdata
fastq_trim_b=${fastq_path}/b_merge
fastq_trim_c=${fastq_path}/c_trim_single
fastq_trim_d=${fastq_path}/d_single_test100

depth='100000'

## fastq ###

# put all fastq files together in ${fastq_raw_a}
# for each sample, merge the R1 and R2 fastq files into a single file
# if you're using single-end reads, you can this step

mkdir -p ${fastq_raw_a}
mv *.gz ${fastq_raw_a}
mkdir -p ${fastq_trim_b}
mkdir -p ${fastq_trim_c}
mkdir -p ${fastq_trim_d} 

while read -r ID;
    do   
#    echo "Processing sample: ${ID}"
    cat ${fastq_raw_a}/${ID}-PCR_R1.fastq.gz ${fastq_raw_a}/${ID}-PCR_R2.fastq.gz > ${fastq_trim_b}/${ID}-PCR_merge.fastq.gz
#    echo "Merged fastq for ${ID}: ${fastq_trim_b}/${ID}-PCR_merge.fastq.gz"
    cd ${fastq_trim_b}
#    echo "Trimming for ${ID}: ${fastq_trim}/${ID}-PCR_merge.fastq.gz"
    $fastp -l 100 -i ${fastq_trim_b}/${ID}-PCR_merge.fastq.gz -o ${fastq_trim_c}/${ID}-PCR_merge_qc.fastq.gz

    done<${SampleList}

for node in $(seq 1 100); do
    echo "Processing node $node..."
    while read -r ID; do
        echo "Trimming R1 for ${ID}: ${fastq_trim_c}/${ID}-PCR_merge_qc.fastq.gz"
        ${seqtk} sample -s${node} ${fastq_trim_c}/${ID}-PCR_merge_qc.fastq.gz ${depth} > ${fastq_trim_d}/${ID}-PCR_merge_qc_s${node}_${depth}.fastq
    done < ${SampleList}
done

