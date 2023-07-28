#!/bin/bash
clear
num_cores=12

rm -rf bcftools
mkdir bcftools

date
echo " "
echo --- Extract Genotypes ---

bcftools mpileup \
  --output-type b \
  --fasta-ref ./reference/GCF_011762595.1_mTurTru1.mat.Y_genomic.fna \
  --annotate FORMAT/DP,FORMAT/AD \
  --max-depth 8000 \
  --skip-indels \
  --threads ${num_cores} \
  ./bam/*.bam | \
bcftools call \
  --variants-only \
  --multiallelic-caller \
  --output-type v \
  --threads ${num_cores} \
  --output bcftools/bcftools.vcf

bcftools stats --verbose \
  bcftools/bcftools.vcf  > \
  ./bcftools/bcftools_vcf_stats.txt

echo " "
date
