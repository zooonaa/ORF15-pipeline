import pandas as pd
import glob
import os
import json


base_dir = os.path.join(os.path.dirname(__file__), "run")  

files = glob.glob("*_VC_final.txt")

for file in files:
    base_name = file.replace("_VC_final.txt", "")

    print(f"processing: {file}")

    try:
        df_raw = pd.read_csv(file, sep="\t", header=None)
        df = df_raw.iloc[:, [1, 2, 4, 5, 12]]
        df.columns = ['Chr', 'Pos', 'Ref', 'Alt', 'Zygosity']
        def classify(z):
            z = str(z).lower()
            if z == 'het':
                return 'het'
            elif z in ['hem', 'hom']:
                return 'hom/hem'
            else:
                return 'other'

        df['Zygosity_class'] = df['Zygosity'].apply(classify)

        summary = df.groupby(['Chr', 'Pos', 'Ref', 'Alt', 'Zygosity_class']).size().unstack(fill_value=0)
        
        for col in ['het', 'hom/hem']:
            if col not in summary.columns:
                summary[col] = 0

        summary['total'] = summary['het'] + summary['hom/hem']
        summary = summary.reset_index()[['Chr', 'Pos', 'Ref', 'Alt', 'het', 'hom/hem', 'total']]

        if summary.empty:
            summary = pd.DataFrame(columns=['Chr', 'Pos', 'Ref', 'Alt', 'het', 'hom/hem', 'total'])

        output_file = f"{base_name}_variant.txt"
        summary.to_csv(output_file, sep="\t", index=False)

        print(f"output: {output_file}")

    except Exception as e:
        print(f"FAIL {file} {e}")
