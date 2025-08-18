import os
import glob
import pandas as pd

def summarize_variant_table(file_pattern: str, output_file: str):
    files = glob.glob(file_pattern)
    all_data = []
    all_samples = []

    print(f"\n🔍 收集檔案：{file_pattern}")
    if not files:
        print("⚠️ 找不到任何 _variant.txt 檔案")
        return

    for file in files:
        sample = file.split("_")[0]
        all_samples.append(sample)
        print(f"  📄 處理 sample: {sample}")

        try:
            df = pd.read_csv(file, sep="\t")

            df['Variant'] = (
                df['Chr'].astype(str) + ":" +
                df['Pos'].astype(str) + ":" +
                df['Ref'].astype(str) + ":" +
                df['Alt'].astype(str)
            )

            # 根據 score 篩選
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
            print(f"    ❌ 錯誤：{e}")

    if not all_data:
        print("⚠️ 所有 sample 都沒有有效變異，未產出統整表")
        return

    merged = pd.concat(all_data)

    matrix = merged.pivot(index='Variant', columns='Sample', values='Result')

    for sample in all_samples:
        if sample not in matrix.columns:
            matrix[sample] = None

    matrix = matrix.dropna(how='all')

    matrix = matrix.reset_index()
    variant_split = matrix['Variant'].str.split(":", expand=True)
    variant_split.columns = ['Chr', 'Pos', 'Ref', 'Alt']
    matrix = pd.concat([variant_split, matrix.drop(columns='Variant')], axis=1)

    sample_cols = sorted([col for col in matrix.columns if col not in ['Chr', 'Pos', 'Ref', 'Alt']])
    matrix = matrix[['Chr', 'Pos', 'Ref', 'Alt'] + sample_cols]

    matrix.to_excel(output_file, index=False)
    print(f"✅ 已輸出：{output_file}")

# 執行範例
summarize_variant_table("*_haplotypecaller_all_variant.txt", "HC_summary_matrix_from_variant_txt.xlsx")
summarize_variant_table("*_mutect2_all_variant.txt", "M2_summary_matrix_from_variant_txt.xlsx")
