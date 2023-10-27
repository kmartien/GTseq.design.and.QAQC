# This script uses samtools to calculate read depth for each individual at each locus,
# then uses the function summarize.depth.at.target.SNPs to calculate summary statistics
# and produce summary graphs.
#
# Requirements:
# 1) you must have samtools installed
# 2) all of the bam files that you want to summarize need to be in the same folder
# 3) you need a bed file that specifies the sites at which you want to calculate read depth

library(tidyverse)
library(dplyr)
library(ggplot2)
library(Mnov.GTseq.data)
library(swfscMisc)
source("R/functions/summarize.depth.at.target.SNPs.R")

project <- "RunMS51"
bam.dir <- paste0("data-raw/bam.files/", project)
data("SWFSC.bed")
loci <- SWFSC.bed
names(loci) <- c("locus", "start", "stop")

####################################
# Open terminal window, navigate to folder with bam files, then paste:

for FILE in *.bam; do samtools depth -b /Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.data/data-raw/Mnov.targetSNP.bed $FILE > ${FILE%.bam}.coverage; done

####################################

cov.files <- list.files(path = bam.dir, pattern = "*.coverage")

depth.sum <- summarize.depth.at.target.SNP(bam.dir, cov.files, loci)

pdf(file = paste0("results-raw/coverage.summary.",project,".pdf"))
depth.sum$plots$mean.cov
depth.sum$plots$med.cov
depth.sum$plots$genos.v.mean.depth
dev.off()  

write.csv(loc.depth, file = paste0("results-raw/", project, ".LOCdepth.at.target.SNPs.csv"))
write.csv(depth.per.ind, file = paste0("results-raw/", project, ".INDdepth.at.target.SNPs.csv"))
save(depth.sum, file = paste0("results-R/", project, ".depth.at.target.SNPs.rda"))
