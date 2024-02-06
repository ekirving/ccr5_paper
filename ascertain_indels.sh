#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2024, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# Ascertain a set of control INDELs with similar characteristics to CCR5delta32 
# rs333 - chr3:46414947-46414978

# CCR5delta32
BP_CCR5=32      # the size of the deletion
AF_CCR5=0.1098  # the EUR frequency

# find INDELS within +/- range
BP_DIFF=5
AF_DIFF=0.02

# normalise chr3 of the 1000 Genomes Project so all multi-allelic sites are merged into a single record
bcftools norm \
  --multiallelics +any \
  -Oz -o 1000g.chr3.norm.vcf.gz \
  data/1000G/vcf/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# get all the biallelic INDELs (and drop the genotypes to make things faster)
bcftools view \
  --types indels \
  --max-alleles 2 \
  --drop-genotypes \
  -Oz -o 1000g.chr3.norm.indels.vcf.gz \
  1000g.chr3.norm.vcf.gz

# count the biallelic INDELs
bcftools view --no-header 1000g.chr3.norm.indels.vcf.gz | wc -l
# 209,120

# get the list of INDELs with a comparable size and frequency
bcftools view \
  --include "abs(strlen(ALT)-strlen(REF))>=(${BP_CCR5}-${BP_DIFF}) & abs(strlen(ALT)-strlen(REF))<=(${BP_CCR5}+${BP_DIFF}) & EUR_AF>=(${AF_CCR5}-${AF_DIFF}) & EUR_AF<=(${AF_CCR5}+${AF_DIFF})" \
  -Oz -o 1000g.chr3.norm.indels.filtered.vcf.gz \
  1000g.chr3.norm.indels.vcf.gz
  
# count the paired INDELs
bcftools view --no-header 1000g.chr3.norm.indels.filtered.vcf.gz | wc -l
# 15

# convert the VCF into PLINK binary format
plink --make-bed --biallelic-only --set-missing-var-ids '@:#' --vcf 1000g.chr3.norm.vcf.gz --out 1000g.chr3.norm

# get a list of the EUR sample IDs
awk '$3=="EUR" {print $1" "$1}' data/1000G/1000G.poplabels > 1000g.eur.samples

# extract the EUR samples
plink --make-bed --bfile 1000g.chr3.norm --keep 1000g.eur.samples --out 1000g.chr3.norm.eur

# now calculate the pair-wise LD within the EUR subpopulation
plink \
  --bfile 1000g.chr3.norm.eur \
  --r2 in-phase \
  --ld-snps `bcftools query --format '%CHROM:%POS,' 1000g.chr3.norm.indels.filtered.vcf.gz` \
  --ld-window-kb 1000 \
  --ld-window 99999 \
  --ld-window-r2 0.8 \
  --out 1000g.chr3.norm.eur.ld

