import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- ClinVar side: all BRCA1 records ---
cv = pd.read_csv("results/clinvar_brca1.tsv", sep="\t", header=None,
                 na_values=["."], keep_default_na=False,
                 names=["pos", "ref", "alt", "sig", "mc", "gene"])
cv = cv.replace("", pd.NA)
cv = cv[cv.gene.fillna("").str.contains(r"(?:^|\|)BRCA1:")].copy()
# First listed molecular consequence, e.g. "SO:0001583|missense_variant"
cv["consequence"] = (cv.mc.fillna("").str.split(",").str[0]
                     .str.split("|").str[-1])

def sig_group(s):
    s = "" if pd.isna(s) else s
    if "onflicting" in s: return "ClinVar Conflicting"
    if "athogenic" in s and "enign" not in s: return "ClinVar Pathogenic / Likely pathogenic"
    if "ncertain" in s: return "ClinVar VUS"
    if "enign" in s: return "ClinVar Benign / Likely benign"
    return "Other"
cv["group"] = cv.sig.apply(sig_group)

def cat(c):
    c = str(c)
    if "frameshift" in c: return "Frameshift"
    if c in ("nonsense", "stop_gained"): return "Nonsense"
    if c in ("splice_donor_variant", "splice_acceptor_variant"): return "Canonical splice site"
    if "missense" in c: return "Missense"
    if "synonymous" in c: return "Synonymous"
    if "inframe" in c: return "In-frame indel"
    if "intron" in c: return "Intronic"
    return "Other (UTR, non-coding, etc.)"
cv["cat"] = cv.consequence.apply(cat)

# --- Population side: observed 1000 Genomes variants only ---
pop = pd.read_csv("results/brca1_master_table.csv")
pop = pop[pop.observed_in_samples]
pop["cat"] = pop.consequence_main.apply(cat)
pop["group"] = "1000 Genomes (observed)"

all_ = pd.concat([cv[["group", "cat"]], pop[["group", "cat"]]], ignore_index=True)
all_ = all_[all_.group != "Other"]
counts = pd.crosstab(all_.group, all_.cat)
print(counts.to_string())
print()
print("BRCA1 ClinVar records in region:", len(cv))
print(cv.group.value_counts())

props = counts.div(counts.sum(axis=1), axis=0) * 100
print()
print(props.round(1).to_string())

ax = props.plot(kind="barh", stacked=True, figsize=(13, 6), colormap="tab20")
ax.set_xlabel("Share of variants (%)")
ax.set_ylabel("")
ax.set_title("Consequence mix: ClinVar classes vs. 1000 Genomes (BRCA1)")
ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
plt.tight_layout()
plt.savefig("figures/fig5_clinvar_vs_population.png", dpi=300)
counts.to_csv("results/consequence_by_group.csv")
