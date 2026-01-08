import csv
import pandas as pd

DATASET_PATH = "datasets/compounds.csv"
FIXED_PATH = "datasets/compounds_fixed.csv"
FILTERED_PATH = "datasets/compounds_filtered.csv"

with open(DATASET_PATH, "r", encoding="utf-8") as f:
    lines = f.readlines()

fixed_rows = []
for i, line in enumerate(lines):
    line = line.strip().replace('""', '"')
    row = line.split(";")
    if i == 0:
        header_len = len(row)
    else:
        if len(row) < header_len:
            row += [""] * (header_len - len(row))
        elif len(row) > header_len:
            row = row[:header_len]
    fixed_rows.append(row)

with open(FIXED_PATH, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f, delimiter=";", quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerows(fixed_rows)

df = pd.read_csv(FIXED_PATH, sep=";", engine="python")
df.columns = df.columns.str.replace('"', "").str.strip()

for col in df.select_dtypes(include="object").columns:
    df[col] = df[col].str.replace('"', "").str.strip()

filtered_df = df[df["Name"].notna() & (df["Name"].str.strip() != "")]
filtered_df.to_csv(FILTERED_PATH, index=False)

print("Columns:", df.columns.tolist())
print("Filtered rows:", len(filtered_df))
print("Filtered dataset saved to:", FILTERED_PATH)
