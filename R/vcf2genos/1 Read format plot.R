library(vcfR)
library(pinfsc50)
library(reshape2)
library(ggplot2)
library(strataG)
library(tidyverse)
library(dplyr)
#library(easyGgplot2)
library(grid)
library(gridExtra)
library(gdata)
library(Hmisc)
library(Mnov.GTseq.data)
source("R/0 Load vcf2genos functions.R")
source("R/functions/geno.smry.by.locus.R")
source("R/functions/remove.rejected.loci.R")
source("R/functions/microhaplot.2.tgt.R")

description = "GTseq.val.all.new"
vcf.or.microhaplot <- "microhaplot"
date = format(Sys.time(), "%Y%b%d")

loc.ann <- read.csv("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.GTseq.data/data-raw/Mnov_gtseq-locus.annotation.csv")

if(vcf.or.microhaplot == "vcf"){
  vcf <- read.vcfR(paste0("vcf/", description, ".vcf"))
  
  ###### (2)VIOLIN AND (3)HEATMAP PLOTS, (4)COUNT BIALLELIC  ##############################
  violin.plot(vcf, paste0("results/", description, "_", date))
  heat.map(vcf, paste0("results/", description, "_", date))
  count.of.biallelic.loci <- sum(is.biallelic(vcf))
  
  ###### (5)CONVERT VCF TO TIDY ###########################################################
  tgt <- create.tgt(vcf)
} else{
  genos <- readRDS(paste0("/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot/", description, ".rds"))[,-1]
  tgt <- microhaplot.2.tgt(genos)
}
###### (6) CREATE DATAFRAME, GET DEPTH BY LOCUS AND #####################################
###### (7) CREATE DATAFRAME WITH COUNTS AND RATIOS, SUMMARIZE GENO TYPES BY LOCUS #######
geno.df <- create.gtype.df(tgt,paste0("results/", description, "_", date))
loc.dpth <- depth.by.locus(tgt,paste0("results/", description, "_", date))
#create.gtype.cts.ratios.df(tgt,paste0("results/", description, "_", date))

###### (8) REMOVE REJECTED LOCI #############
tgt <- remove.rejected.loci(tgt, loc.ann)
mplot_loci_edit <- filter(mplot_loci_edit, status == "Accept") %>% select(locus, minDP, minAR, maxAR, status)

###### (8) PLOT ALLELE COUNTS BY LOCUS ###################################################
plot.allele.ct.by.locus(tgt=tgt, locus.ARs=mplot_loci_edit, fn.prefix = paste0("results/", description,"_", date))

save(tgt, loc.ann, file = paste0("data/", description, "_", date, "_vcf2tgt.rda"))

###### NEXT: USE ALLELE COUNT PLOTS TO ADJUST AR AND minDP VALUES IN #####################
###### "multiplot_loci_edit.csv", THEN MOVE ON TO XXXX.R. ################################

##################### Adjusting min and maxAR ###################################
#################################################################################
#increasing the maxAR brings the blue lines closer to the green line, slope = 1
#increasing the minAR brings the red lines closer to the blue lines
#decreasing them widens the angle

#slope1(upper blue line) is (1-maxAR)/maxAR
#slope2(lower blue line) is maxAR/(1-maxAR)
#slope3(upper red line) is (1-minAR)/minAR
#slope4(lower red line) is minAR/(1-minAR)

# outside of red lines are homozygotes, inside blue lines are heterozygotes,
#in between red and blue are NA's
