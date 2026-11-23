#!/usr/bin/env bash
# MIN85 SNP-site completeness filtering used for the MycoKeys sensitivity analyses.
# Requires vcftools 0.1.16; input was produced by ipyrad v0.9.102 step 7.
set -euo pipefail

BASE="${PPG_DATA_DIR:?Set PPG_DATA_DIR to the PPG analysis directory}"
ASM="MIN85_Xtar_reduced"
SOURCE="${BASE}/${ASM}_outfiles/${ASM}.vcf"
DIR="${BASE}/VCF_Xtar_reduced_ingroup_MAF05"
mkdir -p "${DIR}"

OUT="${DIR}/${ASM}_ingroup_maf05"
SITE="${OUT}_site80"

# Subset to ingroup before calculating MAF; this is essential for matching the manuscript.
vcftools --vcf "${SOURCE}" \
  --remove-indv TW4688_P_cristiferum_Madagascar \
  --remove-indv TW5198_P_reticulatum_USA_IL \
  --maf 0.05 --recode --out "${OUT}"

# Retain SNP sites called in at least 80% of the 104 ingroup samples.
vcftools --vcf "${OUT}.recode.vcf" \
  --max-missing 0.8 --recode --out "${SITE}"

# Site and individual summaries must be separate vcftools calls.
vcftools --vcf "${SITE}.recode.vcf" --missing-site \
  --out "${SITE}_final_site"
vcftools --vcf "${SITE}.recode.vcf" --missing-indv \
  --out "${SITE}_final_indv"

awk -F '\t' '/^#CHROM/ {print "Ingroup samples:", NF-9}' "${SITE}.recode.vcf"
awk '!/^#/ {n++} END {print "Retained SNP records:", n+0}' "${SITE}.recode.vcf"
echo "Expected for the manuscript: 104 ingroup samples and 844 SNP records."
echo "Original post-MAF SNP records: 27,647; MIN85 pre-filter loci: 928."
