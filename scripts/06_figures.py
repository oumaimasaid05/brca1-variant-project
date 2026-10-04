import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", context="talk")
df = pd.read_csv("results/brca1_master_table.csv")
for c in ["ac_1000g", "an_1000g"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")
df["af_calc"] = df.ac_1000g / df.an_1000g

def clinvar_group(s):
    if pd.isna(s): return "Not in ClinVar"
    if "onflicting" in s: return "Conflicting"
    if "Uncertain" in s: return "VUS"
    if "athogenic" in s and "enign" not in s: return "Pathogenic"
    if "enign" in s: return "Benign / Likely benign"
    return "Other"

df["cv_group"] = df.clinvar_sig.apply(clinvar_group)
cv_colors = {"Benign / Likely benign": "#4c9f70", "VUS": "#e0a030",
             "Conflicting": "#9b6bb5", "Pathogenic": "#d1495b",
             "Not in ClinVar": "#b0b0b0", "Other": "#555555"}

# Fig 1: consequence distribution
fig, ax = plt.subplots(figsize=(11, 7))
order = df.consequence_main.value_counts()
sns.barplot(x=order.values, y=order.index, color="#3b6ea5", ax=ax)
ax.set_xscale("log")
ax.set_xlabel("Number of variants (log scale)")
ax.set_ylabel("")
ax.set_title("BRCA1 variant consequences (MANE Select, n=%d)" % len(df))
for i, v in enumerate(order.values):
    ax.text(v * 1.05, i, str(v), va="center", fontsize=12)
plt.tight_layout(); plt.savefig("figures/fig1_consequences.png", dpi=300); plt.close()

# Fig 2: ClinVar classification
fig, ax = plt.subplots(figsize=(9, 6))
cnt = df.cv_group.value_counts()
sns.barplot(x=cnt.index, y=cnt.values, hue=cnt.index,
            palette=cv_colors, legend=False, ax=ax)
ax.set_yscale("log")
ax.set_ylabel("Number of variants (log scale)")
ax.set_xlabel("")
ax.set_title("ClinVar classification of BRCA1 variants")
plt.xticks(rotation=25, ha="right")
for i, v in enumerate(cnt.values):
    ax.text(i, v * 1.1, str(v), ha="center", fontsize=12)
plt.tight_layout(); plt.savefig("figures/fig2_clinvar.png", dpi=300); plt.close()

# Fig 3: allele frequency by consequence (variants with AC > 0)
top = ["intron_variant", "upstream_gene_variant", "synonymous_variant",
       "missense_variant"]
sub = df[(df.consequence_main.isin(top)) & (df.ac_1000g > 0)].copy()
fig, ax = plt.subplots(figsize=(10, 6))
sns.boxplot(data=sub, x="consequence_main", y="af_calc", order=top,
            color="#7fa8d1", ax=ax)
ax.set_yscale("log")
ax.set_ylabel("Allele frequency (1000 Genomes, log scale)")
ax.set_xlabel("")
ax.set_title("Allele frequency by consequence")
plt.xticks(rotation=20, ha="right")
plt.tight_layout(); plt.savefig("figures/fig3_af_by_consequence.png", dpi=300); plt.close()

# Fig 4: lollipop plot of coding variants along the protein
cod = df[df.consequence_main.isin(
    ["missense_variant", "synonymous_variant", "stop_gained"])].copy()
cod["aa_pos"] = pd.to_numeric(
    cod.Protein_position.astype(str).str.extract(r"(\d+)")[0], errors="coerce")
cod = cod.dropna(subset=["aa_pos"])

L = 1863  # BRCA1 protein length (UniProt P38398)
domains = [("RING", 24, 65, "#8ecae6"),
           ("BRCT 1", 1646, 1736, "#ffb703"),
           ("BRCT 2", 1756, 1855, "#fb8500")]  # verify boundaries in UniProt P38398

fig, ax = plt.subplots(figsize=(14, 5))
ax.add_patch(plt.Rectangle((0, -0.35), L, 0.25, color="#dddddd"))
for name, a, b, col in domains:
    ax.add_patch(plt.Rectangle((a, -0.35), b - a, 0.25, color=col))
    ax.text((a + b) / 2, -0.62, name, ha="center", fontsize=11)
for g, sub_g in cod.groupby("cv_group"):
    ax.vlines(sub_g.aa_pos, 0, 1, color=cv_colors[g], alpha=0.6, linewidth=1)
    ax.scatter(sub_g.aa_pos, [1] * len(sub_g), color=cv_colors[g], s=40,
               label=f"{g} (n={len(sub_g)})", zorder=3)
ax.set_xlim(0, L); ax.set_ylim(-0.8, 1.3)
ax.set_yticks([]); ax.set_xlabel("Amino acid position (BRCA1, NP_009225.1)")
ax.set_title("Coding variants along the BRCA1 protein")
ax.legend(loc="upper center", ncol=3, fontsize=10, frameon=False)
plt.tight_layout(); plt.savefig("figures/fig4_lollipop.png", dpi=300); plt.close()

print("Coding variants plotted:", len(cod))
print("Done: figures/fig1 .. fig4")
