# BRCA1 variant analysis and clinical annotation

Personal bioinformatics project, October 2026.

## Question
Which variants occur in BRCA1 in a general-population reference panel, what are their predicted consequences, and which have established clinical significance in ClinVar?

## Data and tools
- 1000 Genomes, release 20190312 (biallelic SNVs and indels, GRCh38, 2,548 samples). Chromosome 17, region 17:43,044,295-43,170,245, sites only
- ClinVar GRCh38 VCF, downloaded 2026-10-04
- Ensembl VEP web tool, archive site jun2026.archive.ensembl.org (GRCh38; MANE Select, HGVS, co-located variants, gnomAD, SIFT, PolyPhen)
- bcftools 1.24; Python 3.11 (pandas, matplotlib, seaborn)
- UniProt P38398 (entry version 285) for domain coordinates: RING 24-65, BRCT 1 1642-1736, BRCT 2 1756-1855

## Workflow
1. `scripts/01_extract_and_filter.sh`: extract the region, keep biallelic records, annotate with ClinVar (exact allele match)
2. Manual step: Ensembl VEP web tool on `data/brca1_filtered.vcf`; TXT output saved as `data/vep_results.txt`
3. `scripts/02_parse_vep.py`, `03_check_dropped.py`, `05_check_unmatched.py`: parsing and data-integrity checks
4. `scripts/04_build_master_table.py`: one row per variant on the BRCA1 MANE Select transcript (NM_007294.4), joined with ClinVar and allele counts
5. Manual step: ClinVar records in the region exported with `bcftools query -r 17:43044295-43170245 -f '%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/MC\t%INFO/GENEINFO\n' data/clinvar.vcf.gz > results/clinvar_brca1.tsv` (the result is committed)
6. `scripts/06_figures.py`, `07_clinvar_vs_population.py`: figures
7. `scripts/08_flag_observed.py`, `09_notable_variants.py`: flag variants carried by no sample; list notable variants

The `data/` folder is not tracked because of file size.

## Results

| Step | Variants |
|---|---|
| Raw region | 3,143 (2,896 SNVs, 247 indels); all PASS and biallelic |
| Outside the BRCA1 MANE Select annotation (excluded) | 1,002 (903 in NBR2, NBR1 and the BRCA1P1 pseudogene; 99 with no gene annotation) |
| On the BRCA1 MANE Select transcript | 2,141 |
| of which carried by at least one sample (AC > 0) | 2,101 |
| of which listed in the callset with AC = 0 | 40 |

Consequences (2,141 variants): intronic 1,728; upstream 204; missense 99; synonymous 48; 3' UTR 33; 5' UTR 7; splice region or polypyrimidine tract 16; stop gained 2; canonical splice site 2; other splice-donor terms 2.

ClinVar classification (2,141 variants): not in ClinVar 1,465; benign 503; likely benign 82; benign/likely benign 22; conflicting 48; uncertain significance 18; pathogenic 3.

Pathogenic variants: c.4357+1G>A (splice donor, 1 carrier), c.34C>T p.Gln12Ter (1 carrier), c.3403C>T p.Gln1135Ter (listed in the callset, no carriers).

Consequence mix by ClinVar class, compared with the 2,101 observed 1000 Genomes variants:

| Group | n | Loss of function | Missense | Intronic |
|---|---|---|---|---|
| ClinVar pathogenic / likely pathogenic | 4,116 | 91.5% | 5.7% | 1.7% |
| ClinVar uncertain significance | 2,491 | 4.0% | 77.7% | 13.2% |
| ClinVar conflicting | 2,924 | 1.0% | 91.1% | 6.8% |
| ClinVar benign / likely benign | 3,732 | 0.3% | 13.1% | 76.1% |
| 1000 Genomes, observed | 2,101 | 0.1% | 4.7% | 80.6% |

Loss of function = frameshift, nonsense and canonical splice-site variants.

## Interpretation
- ClinVar pathogenic BRCA1 variants are overwhelmingly loss-of-function. Only 3 of 2,101 observed population variants are predicted loss-of-function, which is consistent with selection against truncating BRCA1 variants.
- Predicted impact is not clinical classification: one canonical splice-acceptor variant (c.-19-1G>A) is HIGH impact in VEP but classified as uncertain significance in ClinVar.
- Missense variants carry most of the uncertainty: 42 of 98 observed missense variants are uncertain or conflicting in ClinVar. SIFT and PolyPhen disagree for two of the four missense VUS, so in-silico predictions alone are not decisive.
- Many VUS are non-coding and have weak review status (single submitter or no assertion criteria).

## Limitations
- 1000 Genomes is a reference panel, not a patient cohort. Carrying a ClinVar-pathogenic allele here says nothing about an individual's clinical status.
- ClinVar is biased toward variants submitted from clinical testing, so its consequence mix is not a population frequency.
- Only the first listed consequence is used per ClinVar record.
- Allele frequencies are computed as AC/AN because the AF field in the file is rounded.
- 40 variants are listed with AC = 0. Figures 1, 2 and 4 include them; figures 3 and 5 use observed variants only.
- Variants are annotated with existing ClinVar classifications. No ACMG/AMP classification was performed.

## Figures
`figures/fig1_consequences.png`, `fig2_clinvar.png`, `fig3_af_by_consequence.png`, `fig4_lollipop.png`, `fig5_clinvar_vs_population.png`
