import pandas as pd
import glob
import os
import json

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.join(script_dir, "run")
os.makedirs(base_dir, exist_ok=True)

files = glob.glob(os.path.join(base_dir, "*_VC_final.txt"))

if not files:
    print(f"[WARN] No *_VC_final.txt found under: {base_dir}")

for file in files:
    base_name = os.path.basename(file).replace("_VC_final.txt", "")
    print(f"processing: {file}")

    try:
        df_raw = pd.read_csv(file, sep="\t", header=None)

        needed_idx = [1, 2, 4, 5, 12]
        if df_raw.shape[1] <= max(needed_idx):
            raise ValueError(f"Not enough columns in {file}. Found {df_raw.shape[1]} columns.")

        df = df_raw.iloc[:, needed_idx]
        df.columns = ['Chr', 'Pos', 'Ref', 'Alt', 'Zygosity']

        def classify(z):
            z = str(z).strip().lower()
            if z == 'het':
                return 'het'
            elif z in ['hem', 'hom']:
                return 'hom/hem'
            else:
                return 'other'

        df['Zygosity_class'] = df['Zygosity'].apply(classify)

        summary = (
            df.groupby(['Chr', 'Pos', 'Ref', 'Alt', 'Zygosity_class'])
              .size().unstack(fill_value=0)
        )

        for col in ['het', 'hom/hem']:
            if col not in summary.columns:
                summary[col] = 0

        summary['total'] = summary['het'] + summary['hom/hem']
        summary = summary.reset_index()[['Chr', 'Pos', 'Ref', 'Alt', 'het', 'hom/hem', 'total']]

        if summary.empty:
            summary = pd.DataFrame(columns=['Chr', 'Pos', 'Ref', 'Alt', 'het', 'hom/hem', 'total'])

        output_file = os.path.join(base_dir, f"{base_name}_variant.txt")
        summary.to_csv(output_file, sep="\t", index=False)
        print(f"output: {output_file}")

    except Exception as e:
        print(f"FAIL {file} {e}")

