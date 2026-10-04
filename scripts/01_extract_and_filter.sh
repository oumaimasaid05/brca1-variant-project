#!/usr/bin/env bash
# Extract the BRCA1 region from 1000 Genomes (GRCh38), filter, and annotate with ClinVar.
# Run from the project root. Contig names in these files have no "chr" prefix.
set -euo pipefail

URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"

mkdir -p data results figures

# Sites-only extraction of the BRCA1 region (GRCh38)
bcftools view -r 17:43044295-43170245 -G -Oz -o data/brca1_1000G.vcf.gz "$URL"
tabix -p vcf data/brca1_1000G.vcf.gz
rm -f ALL.chr17.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz.tbi

# Keep biallelic records only
bcftools view -m2 -M2 data/brca1_1000G.vcf.gz -Oz -o data/brca1_filtered.vcf.gz
tabix -p vcf data/brca1_filtered.vcf.gz

echo "Raw: $(bcftools view -H data/brca1_1000G.vcf.gz | wc -l)"
echo "Filtered: $(bcftools view -H data/brca1_filtered.vcf.gz | wc -l)"
echo "Rare AF<0.01: $(bcftools view -H -i 'INFO/AF<0.01' data/brca1_filtered.vcf.gz | wc -l)"
echo "Common AF>=0.01: $(bcftools view -H -i 'INFO/AF>=0.01' data/brca1_filtered.vcf.gz | wc -l)"

# Uncompressed VCF for the Ensembl VEP web tool
bcftools view data/brca1_filtered.vcf.gz -Ov -o data/brca1_filtered.vcf

# ClinVar (GRCh38) and exact-allele annotation
wget -nc -P data https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget -nc -P data https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
bcftools annotate -a data/clinvar.vcf.gz \
  -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN \
  --pair-logic exact \
  data/brca1_filtered.vcf.gz -Oz -o data/brca1_clinvar.vcf.gz
tabix -p vcf data/brca1_clinvar.vcf.gz
