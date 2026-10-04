import pandas as pd

df = pd.read_csv("results/brca1_master_table.csv")
df["ac_1000g"] = pd.to_numeric(df["ac_1000g"], errors="coerce")
df["observed_in_samples"] = df["ac_1000g"] > 0
df.to_csv("results/brca1_master_table.csv", index=False)

print("BRCA1 analysis set:", len(df))
print("Observed in >=1 sample (AC>0):", int(df.observed_in_samples.sum()))
print("Listed but not observed (AC=0):", int((~df.observed_in_samples).sum()))
print()
print("Consequences of unobserved variants:")
print(df[~df.observed_in_samples].consequence_main.value_counts())
print()
print("ClinVar classes of unobserved variants:")
print(df[~df.observed_in_samples].clinvar_sig.value_counts(dropna=False))
