import pandas as pd

df = pd.read_csv("results/brca1_master_table.csv")
df["af_calc"] = pd.to_numeric(df.ac_1000g, errors="coerce") / pd.to_numeric(df.an_1000g, errors="coerce")
s = df.clinvar_sig.fillna("")

notable = df[(df.IMPACT == "HIGH")
             | s.str.contains("athogenic") & ~s.str.contains("onflicting")
             | (df.observed_in_samples & s.str.contains("Uncertain"))].copy()

cols = ["pos", "rsid", "HGVSc", "HGVSp", "consequence_main", "IMPACT",
        "ac_1000g", "af_calc", "observed_in_samples",
        "clinvar_sig", "clinvar_revstat", "SIFT", "PolyPhen"]
notable = notable[cols].sort_values("pos")
notable["HGVSc"] = notable.HGVSc.str.split(":").str[-1]
notable["HGVSp"] = notable.HGVSp.str.split(":").str[-1]
notable.to_csv("results/notable_variants.csv", index=False)

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 30)
print(len(notable), "notable variants")
print(notable.to_string(index=False))

print()
obs = df[df.observed_in_samples]
print("Observed missense variants:", (obs.consequence_main == "missense_variant").sum())
print("  of which ClinVar-classified:", ((obs.consequence_main == "missense_variant") & obs.clinvar_sig.notna()).sum())
print("  of which VUS or conflicting:", ((obs.consequence_main == "missense_variant")
      & obs.clinvar_sig.fillna("").str.contains("ncertain|onflicting")).sum())
