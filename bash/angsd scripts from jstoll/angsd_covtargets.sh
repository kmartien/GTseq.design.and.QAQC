#!/bin/bash
#run ANGSD for all sorted and filtered bam files for SNP discovery and determining allele frequencies
#need to run samtools to index all input bam files before running this script
#comparing to loci only from Rapture targets


#SBATCH --time 36:00:00
#SBATCH -p cpu-long
#SBATCH -J "covrestrict_snp_disco"
#SBATCH -o log_%j%x
#SBATCH -c 20
#SBATCH --mem=8000



module load angsd/0.935

nInd=183

minInd=100

outdir=/nese/meclab/Jamie/RAD_GTHI/cm_RAD_2/6_SNPdisco

angsd -bam $outdir/all.bams.txt -out $outdir/covrestrict_0.05_IND100 -P 20 -only_proper_pairs 1 -uniqueOnly 1 -minQ 20 -minMapQ 10 -GL 1 -doMajorMinor 1 -doMaf 2 -minInd $minInd -minMaf 0.05 -rf ./filtered_regions_ANGSD.txt
