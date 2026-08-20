#!/usr/bin/env bash
set -e

IN=${1:-bi_10405.vcf.gz}
KEEP=${2:-config/bi_keep_samples_32.txt}
OUT=${3:-subset_exHybrid.nex}

bcftools view -S "$KEEP" -m2 -M2 -v snps -Oz \
  -o bi_subset.vcf.gz "$IN"

bcftools query -l bi_subset.vcf.gz > bi_samples.txt
bcftools query -f '%REF\t%ALT[\t%GT]\n' bi_subset.vcf.gz > bi_geno.tsv

python3 scripts/vcf_to_nexus.py bi_samples.txt bi_geno.tsv "$OUT"
cat scripts/09_mrbayes_block.txt >> "$OUT"
