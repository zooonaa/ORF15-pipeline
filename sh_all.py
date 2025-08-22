import os, sys
import json

#=================================basic========================================#

script_dir = os.path.dirname(os.path.abspath(__file__))
os.system('mkdir '+'run/')
os.system('cd run/')
path1 = os.path.join(script_dir, "run/")
SampleList = os.path.join(script_dir, "list.txt")
fastq_path = os.path.join(script_dir, "list.txt")


fastq_trim_d=fastq_path+"d_single_test100"

data="_PCR"
tool = '_all'
align = '_bwamem'
threading ='16'
#=======================software direction=====================================#
with open(os.path.join(script_dir, "config.json")) as f:
    config = json.load(f)

bwa = os.path.join(config["bwa"])
picard = os.path.join(config["picard"])
GATK4 = os.path.join(config["GATK4"])
samtools = os.path.join(config["samtools"])
ref_ucsc = os.path.join(config["ref_ucsc"])
ref_fai = os.path.join(config["ref_fai"])
samclip = os.path.join(config["samclip"])

#sample_list = os.path.join(script_dir, config["SampleList"])
#fastq_dir = os.path.join(script_dir, config["fastq_dir"])
#output_dir = os.path.join(script_dir, config["output_dir"])

jar = 'java  -Xmx80g -jar'; #40g/80g (V4_VERSION=80)
python3="/bin/python3.6"
#=============================================================#


with open(SampleList, 'r') as file:
    S_List = file.readlines()
S_List = [i.strip() for i in S_List]
node_range = range(1, 101)
for i in S_List:
    sample_dir = os.path.join(path1, i)
    os.makedirs(sample_dir, exist_ok=True)
    for node in node_range:
        status = f"s{node}"
        status_dir = os.path.join(sample_dir, status)
        os.makedirs(status_dir, exist_ok=True)
        fastq_r1 = f"/{fastq_trim_d}/{i}-PCR_merge_qc_{status}_100000.fastq"
        os.system('mkdir'+' '+i)
        file1 = open(i+'_'+status+'.sh', 'w', encoding='utf-8')
#========================================================#

        file1.write('#!/usr/bin/sh'+'\n')
        file1.write('#SBATCH -A MST109178'+'\n')        # Account name/project number				（科技部計畫編號！，要去國網上面查詢）
        file1.write('#SBATCH -J '+i+'_RPGR'+'\n')         # Job name		 				（可以自己隨便改）
        file1.write('#SBATCH -p ngs186G'+'\n')           # Partition Name 等同PBS裡面的 -q Queue name  			（除非WES,WGS才會去改，寫法可以去國網上PARTITION,,看）
        file1.write('#SBATCH -c 28'+'\n')               # 使用的core數 請參考Queue資源設定 				（跟partition一起）
        file1.write('#SBATCH --mem=186g'+'\n')          # 使用的記憶體量 請參考Queue資源設定				（跟partition一起）
        file1.write('#SBATCH -o '+i+'_'+status+'_'+'out.log'+'\n')          # Path to the standard output file 				（也可以絕對路徑/out.log）
        file1.write('#SBATCH -e '+i+'_'+status+'_'+'err.log'+'\n')         # Path to the standard error ouput file			（同上）
        file1.write('#SBATCH --mail-user=zona890721@gmail.com'+'\n')    # email					（改成自己的email）
        file1.write('#SBATCH --mail-type=FAIL,END'+'\n')  
    

    ### bwa_mem ###
        file1.write(bwa+' mem -t '+threading+' -T 100 -R '+'\'@RG\\tID:'+i+data+align+'\\tLB:'+i+data+align+'\\tSM:'+i+data+align+'\\tPL:ILLUMINA\\\' '+
                        ref_ucsc+' '+fastq_r1+' > '+path1+i+'/'+status+'/'+i+data+align+'.sam'+' '+'&&'+'\n' )
        file1.write('\n' )


    ### Samclip ###

        file1.write(samclip+' --ref '+ref_fai+' < '+path1+i+'/'+status+'/'+i+data+align+'.sam'+' > '+path1+i+'/'+status+'/'+i+data+align+'_samclip_.sam'+' &&'+'\n')
        file1.write('\n' )


    ### picard_SortSam ### 

        file1.write(jar+' '+picard+' '+'SortSam INPUT='+path1+i+'/'+status+'/'+i+data+align+'_samclip_.sam'+' '+'OUTPUT='+path1+i+'/'+status+'/'+i+data+align+'.bam'+
                      ' '+'SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true'+' '+'&&'+'\n') 
        file1.write('\n' )


    ### samtools!


        file1.write(samtools+' view -h -q 60 -b '+path1+i+'/'+status+'/'+i+data+align+'.bam'+' > '+path1+i+'/'+status+'/'+i+data+align+'_fq.bam'+' '+'&&'+'\n')
        file1.write('\n' )

        file1.write(samtools+' index '+path1+i+'/'+status+'/'+i+data+align+'_fq.bam'+' '+'&&'+'\n')
        file1.write('\n' )


    ### picard_MarkDuplicates ###
        file1.write(jar+' '+picard+' '+'MarkDuplicates'+' '+'INPUT='+path1+i+'/'+status+'/'+i+data+align+'_fq.bam'+' '+'OUTPUT='+path1+i+'/'+status+'/'+i+data+align+'.marked.bam'+' '+'METRICS_FILE='+path1+i+'/'+status+'/'+i+data+align+'_metrics'+' '+'VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true'+' '+'&&'+'\n')
        file1.write('\n' )


    ### sh09_picard_SortSam ###
        file1.write(jar+' '+picard+' '+'SortSam'+' '+'INPUT='+path1+i+'/'+status+'/'+i+data+align+'.marked.bam'+' '+'OUTPUT='+path1+i+'/'+status+'/'+i+data+align+'.marked.indexed.bam'+' '+'SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true'+' '+'&&'+'\n')
        file1.write('\n' )
        

    ### GATK_HaplotypeCaller ###
        file1.write(jar+' '+GATK4+' '+'HaplotypeCaller --minimum-mapping-quality 60 -L X:38145300-38145800 -R'+' '+ref_ucsc+' '+'-I'+' '+path1+i+'/'+status+'/'+i+data+align+'.marked.indexed.bam'+' '+'-O'+' '+path1+i+'/'+status+'/'+i+data+align+'.haplotype_region.SnpIndel.vcf.gz'+' '+'-bamout '+path1+i+'/'+status+'/'+i+data+align+'.haplotype_region.bamout.bam'+'\n')
        file1.write('\n' )



    ### GATK_Mutect2 ###
        file1.write(jar+' '+GATK4+' '+'Mutect2 --minimum-mapping-quality 60 -L X:38145300-38145800 -R'+' '+ref_ucsc+' '+'-I'+' '+path1+i+'/'+status+'/'+i+data+align+'.marked.indexed.bam'+' '+'-O'+' '+path1+i+'/'+status+'/'+i+data+align+'.Mutect2_region.vcf.gz'+' '+'-bamout '+path1+i+'/'+status+'/'+i+data+align+'.Mutect2_region.bamout.bam'+'\n')
        file1.write('\n' )

        file1.write('rm '+path1+i+'/'+status+'/'+i+data+align+'.sam'+'\n')
        file1.write('rm '+path1+i+'/'+status+'/'+i+data+align+'_samclip_.sam'+'\n')
        file1.write('rm '+path1+i+'/'+status+'/'+i+data+align+'.bam'+'\n')
        file1.write('rm '+path1+i+'/'+status+'/'+i+data+align+'_fq.bam'+'\n')

        file1.close()
#        os.system('sbatch '+i+'_'+status+'.sh')
    
