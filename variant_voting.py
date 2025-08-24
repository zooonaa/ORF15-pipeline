#!/usr/bin/env python3
import os
import glob
import shutil
import pandas as pd

script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.join(script_dir, "run")
os.makedirs(base_dir, exist_ok=True)

def summarize_variant_table(file_pattern: str, output_basename: str):

    pattern = os.path.join(base_dir, file_pattern)
    files = sorted(glob.glob(pattern))
    all_data = []
    all_samples = []

    print(f"\n[Collect] {pattern}")
    if not files:
        print("  No matching files found.")
        return

    for file in files:
        fname = os.path.basename(file)

        # Robust sample name parsing
        if fname.endswith("_haplotypecaller_all_variant.txt"):
            sample = fname.replace("_haplotypecaller_all_variant.txt", "")
        elif fname.endswith("_mutect2_all_variant.txt"):
            sample = fname.replace("_mutect2_all_variant.txt", "")
        else:
            sample = fname.rsplit("_variant.txt", 1)[0]

        all_samples.append(sample)
        print(f"  Processing sample: {sample}")

        try:
            df = pd.read_csv(file, sep="\t")

            # Required columns
            need_cols = {'Chr', 'Pos', 'Ref', 'Alt', 'het', 'hom/hem'}
            if not need_cols.issubset(df.columns):
                missing = need_cols - set(df.columns)
                raise ValueError(f"missing columns: {sorted(missing)}")

            df['Variant'] = (
                df['Chr'].astype(str) + ":" +
                df['Pos'].astype(str) + ":" +
                df['Ref'].astype(str) + ":" +
                df['Alt'].astype(str)
            )

            df['het'] = pd.to_numeric(df['het'], errors='coerce')
            df['hom/hem'] = pd.to_numeric(df['hom/hem'], errors='coerce')
            df = df.fillna(0)

            df = df[(df['het'] >= 50) | (df['hom/hem'] >= 50)]

            def pick_call(row):
                if row['het'] >= 50:
                    return f"het_{int(row['het'])}"
                elif row['hom/hem'] >= 50:
                    return f"hom/hem_{int(row['hom/hem'])}"
                else:
                    return ""

            df['Result'] = df.apply(pick_call, axis=1)
            df['Sample'] = sample
            df = df[df['Result'] != ""]

            if not df.empty:
                all_data.append(df[['Variant', 'Sample', 'Result']])

        except Exception as e:
            print(f"    FAIL {fname}: {e}")

    if not all_data:
        print("  No valid variants passed the threshold.")
        return

    merged = pd.concat(all_data, ignore_index=True)

    matrix = merged.pivot(index='Variant', columns='Sample', values='Result')

    for s in all_samples:
        if s not in matrix.columns:
            matrix[s] = None

    matrix = matrix.dropna(how='all')

    matrix = matrix.reset_index()
    variant_split = matrix['Variant'].str.split(":", expand=True)
    variant_split.columns = ['Chr', 'Pos', 'Ref', 'Alt']
    matrix = pd.concat([variant_split, matrix.drop(columns='Variant')], axis=1)

    sample_cols = sorted([c for c in matrix.columns if c not in ['Chr', 'Pos', 'Ref', 'Alt']])
    matrix = matrix[['Chr', 'Pos', 'Ref', 'Alt'] + sample_cols]

    out_xlsx = os.path.join(base_dir, output_basename)
    final_output = None
    try:
        import openpyxl  # noqa: F401
        matrix.to_excel(out_xlsx, index=False)
        final_output = out_xlsx
        print(f"  Output (xlsx): {out_xlsx}")
    except Exception:
        out_tsv = os.path.splitext(out_xlsx)[0] + ".tsv"
        matrix.to_csv(out_tsv, sep="\t", index=False)
        final_output = out_tsv
        print(f"  Output (tsv fallback): {out_tsv}")

    if final_output:
        dest = os.path.join(script_dir, os.path.basename(final_output))
        shutil.copy(final_output, dest)
        print(f"  Copied to: {dest}")

def main():
    summarize_variant_table("*_haplotypecaller_all_variant.txt",
                            "HC_summary_matrix_from_variant_txt.xlsx")
    summarize_variant_table("*_mutect2_all_variant.txt",
                            "M2_summary_matrix_from_variant_txt.xlsx")

if __name__ == "__main__":
    main()

