import pandas as pd

path = "data/vep_results.txt"
with open(path) as f:
    n = sum(1 for l in f if l.startswith("##"))
df = pd.read_csv(path, sep="\t", skiprows=n, na_values=["-", "."], low_memory=False)
df["variant"] = df["Location"] + "_" + df["UPLOADED_ALLELE"]

kept = set(df[(df.SYMBOL == "BRCA1") & (df.MANE == "MANE_Select")].variant)
dropped = df[~df.variant.isin(kept)].drop_duplicates("variant").copy()
dropped["pos"] = dropped["Location"].str.split(":").str[1].str.split("-").str[0].astype(int)

print("Dropped variants:", len(dropped))
print(dropped["pos"].describe())
print(dropped["SYMBOL"].value_counts(dropna=False).head(8))
