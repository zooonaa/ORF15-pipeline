import os
import gzip
import time

# 把全部批次的結果merge在一起。


script_dir = os.path.dirname(os.path.abspath(__file__))
path1 = os.path.join(script_dir, "run/")

output_folder = path1 + "output"
os.makedirs(output_folder, exist_ok=True)


log_file_path = os.path.join(path1, "progress.log")

with open(log_file_path, 'w') as log_file:
    log_file.write("Progress log started\n")

for sample_name in os.listdir(path1):
    sample_path = os.path.join(path1, sample_name)
    if os.path.isdir(sample_path):
        start_time = time.time()

        haplotype_output_path = os.path.join(path1, f"{sample_name}_haplotypecaler_all_VC.txt")
        mutect2_output_path = os.path.join(path1, f"{sample_name}_mutect2_all_VC.txt")

        with open(log_file_path, 'a') as log_file:
            log_file.write(f"Processing sample: {sample_name}\n")

        print(f"Processing sample: {sample_name}...")

        with open(haplotype_output_path, 'w') as haplotype_output_file, \
             open(mutect2_output_path, 'w') as mutect2_output_file:
            
            for sub_folder in os.listdir(sample_path):
                sub_folder_path = os.path.join(sample_path, sub_folder)
                if os.path.isdir(sub_folder_path) and sub_folder.startswith("s"):
                    print(f"Processing folder: {sub_folder} in sample: {sample_name}")
                    with open(log_file_path, 'a') as log_file:
                        log_file.write(f"  Processing folder: {sub_folder}\n")

                    haplotype_pattern = f"{sample_name}_PCR_bwamem.haplotype_region.SnpIndel.vcf.gz"
                    mutect2_pattern = f"{sample_name}_PCR_bwamem.Mutect2_region.vcf.gz"

                    haplotype_file_path = os.path.join(sub_folder_path, haplotype_pattern)
                    if os.path.exists(haplotype_file_path):
                        print(f"  Found haplotype file: {haplotype_file_path}")
                        with gzip.open(haplotype_file_path, 'rt') as haplotype_file:
                            lines = haplotype_file.readlines()[111:]  
                            if lines:
                                for line in lines:
                                    haplotype_output_file.write(f"{sub_folder}\t{line}")
                            else:
                                haplotype_output_file.write(f"{sub_folder}\t0\n")
                    else:
                        with open(log_file_path, 'a') as log_file:
                            log_file.write(f"    Haplotype file not found: {haplotype_file_path}\n")
                        print(f"  Haplotype file not found: {haplotype_file_path}")

                    mutect2_file_path = os.path.join(sub_folder_path, mutect2_pattern)
                    if os.path.exists(mutect2_file_path):
                        print(f"  Found Mutect2 file: {mutect2_file_path}")
                        with gzip.open(mutect2_file_path, 'rt') as mutect2_file:
                            lines = mutect2_file.readlines()[129:]  
                            if lines:
                                for line in lines:
                                    mutect2_output_file.write(f"{sub_folder}\t{line}")
                            else:
                                mutect2_output_file.write(f"{sub_folder}\t0\n")
                    else:
                        with open(log_file_path, 'a') as log_file:
                            log_file.write(f"    Mutect2 file not found: {mutect2_file_path}\n")
                        print(f"  Mutect2 file not found: {mutect2_file_path}")

        elapsed_time = time.time() - start_time
        print(f"Finished processing sample: {sample_name} in {elapsed_time:.2f} seconds.")
        with open(log_file_path, 'a') as log_file:
            log_file.write(f"Finished processing sample: {sample_name} in {elapsed_time:.2f} seconds.\n")

print("step1 done!")



def process_percentage_file(input_suffix, output_suffix):
    for filename in os.listdir(path1):
        if filename.endswith(input_suffix):
            patient_id = filename.split(input_suffix)[0]
            input_path = os.path.join(path1, filename)
            output_path = os.path.join(path1, f"{patient_id}{output_suffix}")

            print(f"processing: {filename}")

            output_lines = []
            variant_found = False  

            with open(input_path, "r") as file:
                for line in file:
                    if line.startswith("#"):
                        continue

                    columns = line.strip().split("\t")
                    if len(columns) >= 11:
                        details = columns[10]  
                        try:
                            values = details.split(":")[1].split(",")
                            numerator = int(values[1])  
                            denominator = int(values[0]) + numerator  
                            proportion_str = f"{numerator}/{denominator}"
                            percentage = (numerator / denominator) * 100 

                            if percentage < 20:
                                category = "low"
                            elif 20 <= percentage < 60:
                                category = "het"
                            elif percentage >= 60:
                                category = "hem"
                            else:
                                category = "-"

                            variant_found = True

                            new_columns = columns[:10] + [proportion_str, f"{percentage:.1f}", category]
                            output_lines.append("\t".join(new_columns))
                        except (IndexError, ValueError):
                            new_columns = columns[:10] + ["ERROR", "ERROR", "ERROR"]
                            output_lines.append("\t".join(new_columns))

            if not variant_found:
                output_lines.append("0")  

            with open(output_path, "w") as out_file:
                out_file.write("\n".join(output_lines))

            print(f"DONE: {output_path}")

process_percentage_file("_haplotypecaler_all_VC.txt", "_haplotypecaller_all_VC_final.txt")

process_percentage_file("_mutect2_all_VC.txt", "_mutect2_all_VC_final.txt")

print("step2 done!")

###============================================================================###
""""
##個別細節
from collections import defaultdict

files = [f for f in os.listdir(path1) if f.endswith("_mutect2_all_VC_final.txt")]

for file in files:
    patient_id = file.split("_")[0]  
    input_file = os.path.join(path1, file)
    output_file = os.path.join(output_folder, f"{patient_id}_mutect2_statistics.txt")

    variant_dict = defaultdict(lambda: defaultdict(str))

    with open(input_file, "r") as f:
        for line in f:
            columns = line.strip().split("\t")
            if len(columns) < 8:
                continue  

            batch = columns[0]  
            variant = f"{columns[1]}:{columns[2]}:{columns[4]}:{columns[5]}"  
            status = columns[-1]  

            if status == "low":
                status = "-"

            variant_dict[variant][batch] = status

    ranges = [(f"s{i+1}", f"s{i+10}") for i in range(0, 100, 10)] + [("s1", "s100")]

    with open(output_file, "w") as out_f:
        headers = ["Variant"] + [f"{r[0]}~{r[1]}_{stat}" for r in ranges for stat in ["hem", "het", "-"]]
        out_f.write("\t".join(headers) + "\n")
        
        for variant, batch_data in variant_dict.items():
            row = [variant]
            for r_start, r_end in ranges:
                hem_count = het_count = dash_count = 0
                for i in range(int(r_start[1:]), int(r_end[1:]) + 1):
                    batch = f"s{i}"
                    status = batch_data.get(batch, "-")
                    if status == "hem":
                        hem_count += 1
                    elif status == "het":
                        het_count += 1
                    elif status == "-":
                        dash_count += 1
                row.extend([str(hem_count), str(het_count), str(dash_count)])
            out_f.write("\t".join(row) + "\n")

    print(f"DONE {patient_id} output: {output_file}")

print("step3 M2 done!")



files = [f for f in os.listdir(path1) if f.endswith("_haplotypecaller_all_VC_final.txt")]

for file in files:
    patient_id = file.split("_")[0] 
    input_file = os.path.join(path1, file)
    output_file = os.path.join(output_folder, f"{patient_id}_haplotypecaller_statistics.txt")

    variant_dict = defaultdict(lambda: defaultdict(str))

    with open(input_file, "r") as f:
        for line in f:
            columns = line.strip().split("\t")
            if len(columns) < 8:
                continue  

            batch = columns[0]  
            variant = f"{columns[1]}:{columns[2]}:{columns[4]}:{columns[5]}"  
            status = columns[-1]  

            if status == "low":
                status = "-"

            variant_dict[variant][batch] = status

    ranges = [(f"s{i+1}", f"s{i+10}") for i in range(0, 100, 10)] + [("s1", "s100")]

    with open(output_file, "w") as out_f:
        headers = ["Variant"] + [f"{r[0]}~{r[1]}_{stat}" for r in ranges for stat in ["hem", "het", "-"]]
        out_f.write("\t".join(headers) + "\n")
        
        for variant, batch_data in variant_dict.items():
            row = [variant]
            for r_start, r_end in ranges:
                hem_count = het_count = dash_count = 0
                for i in range(int(r_start[1:]), int(r_end[1:]) + 1):
                    batch = f"s{i}"
                    status = batch_data.get(batch, "-")
                    if status == "hem":
                        hem_count += 1
                    elif status == "het":
                        het_count += 1
                    elif status == "-":
                        dash_count += 1
                row.extend([str(hem_count), str(het_count), str(dash_count)])
            out_f.write("\t".join(row) + "\n")

    print(f"DONE {patient_id} output: {output_file}")

print("steo3 Haplotypecaller done!")
"""