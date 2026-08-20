# Onychodactylus GBS analysis scripts

Scripts for the GBS analyses used in the Korean *Onychodactylus* study.

## Script order

1. `01_stacks.sh` - de novo Stacks assembly and SNP export from cleaned 100-bp FASTQ files.
2. `02_structure_prepare.sh` - SNP subsampling and STRUCTURE input conversion.
3. `03_structure_run.sh` - STRUCTURE analyses for K = 1-10.
4. `04_pca.sh` - missing-data filtering, marker pruning, exclusion of *O. sillanus*, and PCA.
5. `05_pca_plot.R` - PCA figure.
6. `06_fst.R` - pairwise Weir-Cockerham FST, permutation test, and BH correction.
7. `07_cline.R` - nuclear cline, IBD-like models, mtDNA cline, AICc, and Fig. 5.
8. `08_bi_prepare.sh` / `vcf_to_nexus.py` - VCF to NEXUS conversion for Bayesian inference.
9. `09_mrbayes_block.txt` - MrBayes settings for the variable-SNP matrix.

## Marker sets

The analyses used different SNP matrices.

- Stacks `populations.snps.vcf`: 241,668 SNPs in the retained VCF.
- STRUCTURE: 24,163 SNPs.
- PCA: 7,623 SNPs and 76 individuals after excluding *O. sillanus*.
- Main Bayesian tree: 10,405 variable SNP characters and 32 individuals from eight populations.
- FST: the analysis script and output matrices were retained, but the original PLINK dosage input was not.

## Inputs

All public input files use the sample and population codes reported in the manuscript. `config/popmap.tsv` contains the sample-to-population assignments used by the scripts.

For the cline analysis, `structure_q_mt.tsv` must contain `ind`, `pop`, `qOK`, and `mt`; `coordinates.tsv` must contain `pop`, `lat`, and `lon`.

For FST, `gbs_pruned.raw` is a PLINK `--recode A` dosage file.

The exact raw-read preprocessing command used before Stacks was not retained, so raw FASTQ preprocessing is not included in this repository. The Stacks workflow begins from the cleaned 100-bp FASTQ files used in the de novo assembly.

The exact historical STRUCTURE SNP list and the original FST dosage file were not retained. The STRUCTURE preparation script therefore generates a 24,163-SNP subset matching the recorded marker count, but not the exact historical marker identities.

Because the de novo Stacks loci do not have reference-genome coordinates, loci are placed on a pseudochromosome before the PCA marker-pruning step.
