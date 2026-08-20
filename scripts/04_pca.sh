#!/usr/bin/env bash
set -e

VCF=${1:-stacks/populations.snps.vcf.gz}
REMOVE=${2:-config/remove_Osillanus.txt}

bcftools view -h "$VCF" \
  | grep -v '^##contig=' \
  | awk '/^#CHROM/{print "##contig=<ID=1,length=1000000000>"} {print}' \
  > pseudo.header

bcftools view -H "$VCF" \
  | awk 'BEGIN{OFS="\t"}{$1=1;$2=NR;print}' \
  > pseudo.body

cat pseudo.header pseudo.body | bgzip -c > all.pseudo.vcf.gz
rm -f pseudo.header pseudo.body

plink2 --vcf all.pseudo.vcf.gz --double-id --geno 0.2 \
  --make-pgen --out pca_filtered

plink2 --pfile pca_filtered \
  --indep-pairwise 50 5 0.2 --out pca_prune

plink2 --pfile pca_filtered --extract pca_prune.prune.in \
  --make-bed --out pca_pruned

plink --bfile pca_pruned --remove "$REMOVE" \
  --pca 20 --out pca_no_Osillanus
