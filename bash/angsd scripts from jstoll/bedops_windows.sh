#!/bin/bash
#make sure conda env with bedops is active
#generate bed file of SNPd ensity in 300bp sliding windows along the genome

bedops --chop 300 --stagger 280 -x <(awk -vOFS="\t" '{ print $1, $2-1, $3; }' Brazil_chromextent.bed | sort-bed -) | bedmap --echo --count --delim '\t' - <(awk -vOFS="\t" '{ print $1, $2-1, $2; }' brazil_SNP_mafs.bed | sort-bed -) > Brazil_snpdensity.bed
