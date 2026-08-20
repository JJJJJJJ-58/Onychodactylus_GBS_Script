#!/usr/bin/env bash
set -e

VCF=${1:-stacks/populations.snps.vcf}
POPMAP=${2:-config/popmap.tsv}
N_SNP=${N_SNP:-24163}

if [[ "$VCF" == *.gz ]]; then
  INVCF="$VCF"
else
  bgzip -c "$VCF" > populations.snps.vcf.gz
  INVCF=populations.snps.vcf.gz
fi

tabix -f -p vcf "$INVCF"

bcftools query -f '%CHROM\t%POS\n' "$INVCF" \
  | shuf -n "$N_SNP" \
  | sort -t $'\t' -k1,1 -k2,2n > structure.keep.pos

bcftools view -R structure.keep.pos -Oz \
  -o structure_24163.vcf.gz "$INVCF"

plink --vcf structure_24163.vcf.gz \
  --double-id --allow-extra-chr --recode structure --out strc1

awk -v OFS='\t' 'NR==FNR{p[$1]=$2; next}
FNR<=2{next}
{
  id=$2; pop=(id in p ? p[id] : -9);
  printf "%s\t%s", id, pop;
  for(i=7;i<=NF;i++) printf "\t%s", $i;
  print "";
}' "$POPMAP" strc1.recode.strct_in > strc1.for_structure.txt
