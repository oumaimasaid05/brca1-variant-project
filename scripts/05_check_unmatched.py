import io, subprocess
import pandas as pd

fmt = "%CHROM\t%POS\t%REF\t%ALT\n"
out = subprocess.run(["bcftools", "query", "-f", fmt, "data/brca1_clinvar.vcf.gz"],
                     capture_output=True, text=True, check=True).stdout
cv = pd.read_csv(io.StringIO(out), sep="\t", header=None, keep_default_na=False,
                 names=["chrom", "pos", "ref", "alt"])
cv["key"] = cv.chrom.astype(str) + ":" + cv.pos.astype(str) + ":" + cv.ref + ":" + cv.alt

path = "data/vep_results.txt"
with open(path) as f:
    n_meta = sum(1 for line in f if line.startswith("##"))
vep = pd.read_csv(path, sep="\t", skiprows=n_meta, na_values=["-", "."],
                  keep_default_na=False, low_memory=False)
m = vep[(vep.SYMBOL == "BRCA1") & (vep.MANE == "MANE_Select")].copy()

start = m.Location.str.split(":").str[1].str.split("-").str[0].astype(int)
al = m.UPLOADED_ALLELE.str.split("/", expand=True)
ref, alt = al[0], al[1]
trimmed = (ref.str.len() != alt.str.len()) & (ref.str[0] == alt.str[0])
m["pos"] = start - trimmed.astype(int)
m["key"] = "17:" + m.pos.astype(str) + ":" + ref + ":" + alt
m["ref"], m["alt"] = ref, alt

bad = m[~m.key.isin(set(cv.key))].copy()
bad["type"] = ((bad.ref.str.len() == 1) & (bad.alt.str.len() == 1)).map(
    {True: "SNV", False: "indel/other"})

print("Unmatched:", len(bad))
print(bad["type"].value_counts())
print()
print(bad[["Location", "UPLOADED_ALLELE", "Consequence", "pos", "key"]]
      .head(15).to_string())
