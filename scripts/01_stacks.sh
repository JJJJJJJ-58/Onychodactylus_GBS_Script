#!/usr/bin/env bash
set -e

POPMAP=${1:-config/popmap.tsv}

mkdir -p stacks logs

i=1
for fq in *.clean.L100.fq.gz; do
  s=${fq%.clean.L100.fq.gz}
  ustacks -f "$fq" -o stacks -i $i --name "$s" -m 3 -M 2 -p 32 \
    2>&1 | tee "logs/ustacks_${s}.log"
  i=$((i+1))
done

cstacks -P stacks -M "$POPMAP" -n 2 -p 32 | tee logs/cstacks.log
sstacks -P stacks -M "$POPMAP" -p 32 | tee logs/sstacks.log
tsv2bam -P stacks -M "$POPMAP" -t 32 | tee logs/tsv2bam.log
gstacks -P stacks -M "$POPMAP" -t 32 | tee logs/gstacks.log

populations -P stacks -M "$POPMAP" -t 32 \
  -r 0.7 --min-maf 0.05 --write-single-snp --vcf \
  | tee logs/populations.log
