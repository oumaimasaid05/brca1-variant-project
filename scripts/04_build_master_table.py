import io
import subprocess
import pandas as pd

# 1. ClinVar annotations + frequencies exported by bcftools
fmt = ("%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AC\t%INFO/AN\t"
       "%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n")
out = subprocess.run(["bcftools", "query", "-f", fmt, "data/brca1_clinvar.vcf.gz"],
                     capture_output=True, text=True, check=True).stdout
cv = pd.read_csv(io.StringIO(out), sep="\t", header=None, keep_default_na=False,
                 names=["chrom", "pos", "ref", "alt", "af", "ac", "an",
                        "clinvar_sig", "clinvar_revstat", "clinvar_disease"])
for c in ["af", "ac", "an", "clinvar_sig", "clinvar_revstat", "clinvar_disease"]:
    cv[c] = cv[c].replace(".", pd.NA)
cv["key"] = (cv.chrom.astype(str) + ":" + cv.pos.astype(str) + ":"
             + cv.ref + ":" + cv.alt)

# 2. VEP annotations, BRCA1 MANE Select transcript only
path = "data/vep_results.txt"
with open(path) as f:
    n_meta = sum(1 for line in f if line.startswith("##"))
vep = pd.read_csv(path, sep="\t", skiprows=n_meta, na_values=["-", "."],
                  keep_default_na=False, low_memory=False)
m = vep[(vep.SYMBOL == "BRCA1") & (vep.MANE == "MANE_Select")].copy()

# Location start = original VCF position; UPLOADED_ALLELE = original REF/ALT
m["pos"] = m.Location.str.split(":").str[1].str.split("-").str[0].astype(int)
al = m.UPLOADED_ALLELE.str.split("/", expand=True)
m["key"] = "17:" + m.pos.astype(str) + ":" + al[0] + ":" + al[1]

# 3. Join
keep = ["key", "pos", "Consequence", "IMPACT", "HGVSc", "HGVSp",
        "Protein_position", "Amino_acids", "Existing_variation",
        "SIFT", "PolyPhen", "gnomADg_AF"]
master = m[keep].merge(
    cv[["key", "af", "ac", "an", "clinvar_sig", "clinvar_revstat",
        "clinvar_disease"]],
    on="key", how="left", indicator=True)

print("Variants in master table:", len(master))
print("Matched to annotated VCF:", (master["_merge"] == "both").sum())
master = master.drop(columns="_merge")
master = master.rename(columns={"af": "af_1000g", "ac": "ac_1000g",
                                "an": "an_1000g", "gnomADg_AF": "af_gnomad",
                                "Existing_variation": "rsid"})
master["consequence_main"] = master["Consequence"].str.split(",").str[0]
master.to_csv("results/brca1_master_table.csv", index=False)

print()
print(master["clinvar_sig"].value_counts(dropna=False))

print()
print("=== Pathogenic / likely pathogenic ===")
s = master.clinvar_sig.fillna("")
sel = master[s.str.contains("athogenic") & ~s.str.contains("onflicting")]
print(sel[["pos", "consequence_main", "HGVSc", "HGVSp", "ac_1000g", "an_1000g",
           "af_1000g", "clinvar_sig", "clinvar_revstat"]].to_string())
