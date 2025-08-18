import pandas as pd
import glob
import os

# 找到所有符合格式的 txt 檔
files = glob.glob("*_VC_final.txt")

for file in files:
    base_name = file.replace("_VC_final.txt", "")

    print(f"🔄 處理中：{file}")

    try:
        # 沒有欄位名稱 → 手動抓欄位
        df_raw = pd.read_csv(file, sep="\t", header=None)

        # 使用第 2, 3, 5, 6, 13 欄（index = 1,2,4,5,12）
        df = df_raw.iloc[:, [1, 2, 4, 5, 12]]
        df.columns = ['Chr', 'Pos', 'Ref', 'Alt', 'Zygosity']

        # 分類 zygosity
        def classify(z):
            z = str(z).lower()
            if z == 'het':
                return 'het'
            elif z in ['hem', 'hom']:
                return 'hom/hem'
            else:
                return 'other'

        df['Zygosity_class'] = df['Zygosity'].apply(classify)

        # 統計
        summary = df.groupby(['Chr', 'Pos', 'Ref', 'Alt', 'Zygosity_class']).size().unstack(fill_value=0)
        
        # 確保兩個類別都在
        for col in ['het', 'hom/hem']:
            if col not in summary.columns:
                summary[col] = 0

        # 計算 total
        summary['total'] = summary['het'] + summary['hom/hem']
        summary = summary.reset_index()[['Chr', 'Pos', 'Ref', 'Alt', 'het', 'hom/hem', 'total']]

        # 如果處理後沒有任何 variant，也輸出空檔案
        if summary.empty:
            summary = pd.DataFrame(columns=['Chr', 'Pos', 'Ref', 'Alt', 'het', 'hom/hem', 'total'])

        # 輸出
        output_file = f"{base_name}_variant.txt"
        summary.to_csv(output_file, sep="\t", index=False)

        print(f"✅ 完成輸出：{output_file}")

    except Exception as e:
        print(f"❌ 處理 {file} 時出錯：{e}")
