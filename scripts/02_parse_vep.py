import pandas as pd

path = "data/vep_results.txt"

with open(path) as f:
    n_meta = sum(1 for line in f if line.startswith("##"))

df = pd.read_csv(path, sep="\t", skiprows=n_meta,
                 na_values=["-", "."], low_memory=False)
df = df.rename(columns={"#Uploaded_variation": "Uploaded_variation"})
df["variant"] = df["Location"] + "_" + df["UPLOADED_ALLELE"]

print("Total VEP rows (variant x transcript):", len(df))
print("Unique variants:", df["variant"].nunique())

mane = df[(df["SYMBOL"] == "BRCA1") & (df["MANE"] == "MANE_Select")].copy()
print("Variants with a BRCA1 MANE Select row:", mane["variant"].nunique())
print("Duplicated variants:", mane["variant"].duplicated().sum())

mane["consequence_main"] = mane["Consequence"].str.split(",").str[0]
print()
print(mane["consequence_main"].value_counts())
print()
print(mane["IMPACT"].value_counts())
print()
print(mane["CLIN_SIG"].value_counts(dropna=False).head(15))
