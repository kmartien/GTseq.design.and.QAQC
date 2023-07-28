library(reshape2)
library(ggplot2)
library(strataG)
library(tidyverse)
library(dplyr)
source("R/Freebayes.miscalls.r")
source("R/multiplot_pdf.r")

description = "Dcor_gtseq_LocNum_noMNP_nocomplex_101922"
date = format(Sys.time(), "%Y%b%d")

# Note that this assumes that the date on the saved file you're loading is today's
# date. If it was created previously, replace the file name below with the
# correct one.
load(paste0("data/", description, "_", date, "_filtered-tgt.rda"))


###### (12) FIND ND RECALL SPECIFIC GENOTYPES ############################################
##########################################################################################
# Find and re-call or remove specific genotypes (see locus notes)
# Make a data frame where rows are IDs and columns are genotypes with allele counts.
# First, we will select just the columns that we are going to need, and then name the loci
# with the positions in there.

tmp <- tgt %>%
  select(Indiv, locus, gt, depth.x, depth.y, ratio.x, ratio.y) %>%
  unite(gt_counts, depth.x, depth.y, sep = "/") %>%
  unite(ratios, ratio.x, ratio.y, sep = "/") %>%
  unite(gtypes, gt, gt_counts, ratios, sep = " ")
# then use spreading operation to convert to a matrix of samples by loci
wide <- tmp %>%
  spread(key = locus, value = gtypes)
# separate genotypes into 2 columns/locus
#gdata <- cbind(wide[, 1], alleleSplit(wide[, -1], sep = "/"))


gdata <- as.data.frame(wide)
# look at a small part of that
gdata[1:10, 1:10]
# write csv file
write.csv(gdata, paste(description, "_", date, "_recalled_genotypes&ratios.csv", sep = ""),
          row.names = FALSE)


###### (13) REMOVE IND SAMPLE GENOTYPES ###################################################
##########################################################################################

# Remove individual sample genotypes (coded by Amy Van Cise)
# generate new .csv file (genotype_changes_DATE.csv)
# use the _freebayes_miscalls.csv file, insert gt column before gtypes,
# review genotypes in plots, and recall by adding corrected genotype into gt column
# column headers include CHROM, Indiv, gt (for new genotype).
# Replacement genotypes should be in same format (e.g., A/G),
# and excluded genotypes should be blank or NA.

exclude_genotypes <- read.csv("genotype_changes.csv", header = TRUE,
                              stringsAsFactors=FALSE)


# Combine CHROM and Indiv first:
exclude_genotypes <- exclude_genotypes %>% unite(locus_Indiv, locus, Indiv, sep = "_")

tgt <- tgt %>% unite(locus_Indiv, locus, Indiv, sep = "_", remove = FALSE) %>%
  merge(exclude_genotypes, by="locus_Indiv", all.x=TRUE)%>%
  mutate(gt = ifelse(locus_Indiv %in% exclude_genotypes$locus_Indiv,
                     gt.y, gt.x)) %>%
  select(-c(gt.x,gt.y))
# clean up empty's (currently include "NA", "." and ""). Change all to NA
tgt$gt <- ifelse(tgt$gt=="",NA,tgt$gt)
tgt$gt <- ifelse(tgt$gt==".",NA,tgt$gt)
tgt$depth <- tgt$depth.x + tgt$depth.y
# tgt5 should now have the corrected genotypes
# Check to make sure genotypes were changed appropriately:

changed_genotypes <- tgt[which(tgt$locus_Indiv %in% exclude_genotypes$locus_Indiv),
                         c(1,length(tgt))]
write.csv(changed_genotypes, paste(description, "_", date, "_changed_genotypes.csv", sep = ""),
          row.names = FALSE)

#That should give you a table with CHROM_Indiv and the genotype column for only the rows
#that were in exclude_genotypes, so that you can check to make sure that the right rows
#ended up in the right place.

#Check again to make sure that you didn't miss any freebayes miscalls (1-minAR) to 1 ratio, hets)
gdata <- freebayes.miscalls(tgt, paste("results/", description, "_", date, "_freebayes_miscalls_recheck.csv", sep = ""))

###### (13b) Check for Contaminated individuals  ####################################################
##########################################################################################
# Here we are going to make a list of LabID's that have Na/Na, with coverage >20 and
# ratios between 0.15 and 0.85, which is our criteria for potential contamination
# two lists are exported, one with all individuals and loci that meet that criteria
# and the other is a list by lab ID and frequency of meeting that criteria
# there is also a plot to show the histogram of the frequency, and a pdf that plots
# the Frequency of Na/Na by LabID
NaList <- tgt %>%
  select(Indiv, locus, depth.x, depth.y, ratio.x, ratio.y, depth, gt) %>%
  separate(Indiv, c("species", "id"), sep = "or")
NaList <-filter(NaList, gt == "NA/NA")
NaList <-filter(NaList, depth >= 20)
NaList <-filter(NaList, ratio.x > .15 & ratio.x < .85 & ratio.y >.15 & ratio.y < .85)
NaList$species <- NULL

NaList_Freq<- as.data.frame(table(NaList$id))
FreqNaList<- as.data.frame(table(NaList_Freq$Freq))
barplot(FreqNaList$Freq, names.arg = FreqNaList$Var1)
filter(NaList_Freq$Freq > 25)

write.csv(NaList, paste0(description,"_", date, "_NaList.csv"),
          row.names = FALSE)
write.csv(NaList_Freq, paste0(description,"_", date, "_NaList_Freq.csv"),
          row.names = FALSE)


pdf(paste(description, "_NonLowCovNA_Freq","_", date, ".pdf",sep=""), width = 80, height = 8.5)
ggplot(NaList_Freq, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", color = "black", fill = "grey") +
  labs(title = "Frequency of Na/Na(cov >20,ratios >.15 & < .85) by LabID\n", x = "\nLabID", y = "Frequency\n") +
  theme(axis.text.x = element_text(angle=90))
dev.off()

###### (14) TRANSFORM INTO DATAFRAME FOR EXPORT ##########################################
##########################################################################################


# Make data frame (use tgt4 if individual genotype change section (13) was skipped)
tmp <- tgt %>%
  select(Indiv, locus, gt) %>%
  separate(Indiv, c("species", "id"), sep = "or")
tmp$species <- NULL
#This is specific to DCOR data, to delete the "Dcor" at the begining of the sample ID's.
#tmp4 <- tmp %>%
# unite(loc_pos, locus, POS, sep = "_")


# then use spreading operation to convert to a matrix of samples by loci
wide <- tmp %>%
  spread(key = locus, value = gt)
# separate genotypes into 2 columns/locus
gdata <- cbind(wide[, 1], alleleSplit(wide[, -1], sep = "/"))
gdata <- as.data.frame(gdata)
# look at a small part of that
gdata[1:10, 1:10]
# write csv file
write.csv(gdata, paste(description, "_", date, "_final_genotable.csv", sep = ""), row.names = FALSE)


###### (15) REPLOT FINAL DATASET FOR RECORDS #############################################
##########################################################################################

plot.allele.ct.by.locus(tgt, locus_AR_minDP, 
                        fn.prefix = paste0("results/", description, "_", date, "_FINAL_allele_multiplots.pdf"), 
                        smallest.plot.lim = 50)


###### (16) CREATE DATAFRAME FOR DATABASE IMPORT #########################################
##########################################################################################

#convert full dataset to reduced column dataframe for import to database (1 row per genotype).

# Then make data frame of reduced columns, and remove "species" from beginning of sample IDs
# (to leave just LABIDs)
tmp6 <- tgt %>%
  select(Indiv, locus, gt, depth.x, depth.y) %>%
  separate(Indiv, c("species", "id"), sep = "r")

# This is specific to Dcor data, to delete the "Dcor" at the begining of the sample ID's.
tmp6$species <- NULL
write.csv(tmp6, paste0(description,"_", date, "_final_genotype_data_for_DB.csv"),
          row.names = FALSE)
# saving the data this way gets around issue in Access of allowing only ~250 columns.
# Data will need to be converted to a matrix of samples by loci
# (and split columns for alleles) either using a cross-tab table in Access or an R script
# (e.g., #14 above).

# *************** Not sure how the duplicate sample imports will work ********************


###### (17) QA/QC STEPS ##################################################################
##########################################################################################

# Identify duplicates and summarize loci and individual information
# add a strata column to the csv file before importing it
# ***this section of code can be copied to a separate r script and run on ahi

QAQC <- read.csv(file = "MS28_MS29_MS30_MS37_hatchlings_022223_final_genotable.csv")
QAQC_gtype <-df2gtypes(QAQC, ploidy = 2, id.col = 1, strata.col = 2, loc.col = 3,
                       description = (paste(description,"_",date,"_qaqc_check")))

smry <- summarizeLoci(QAQC_gtype)
write.csv(smry, paste(description, "_", date,"_locus_summary.csv", sep = "",
                      row.names = TRUE))


smry_ind <-summarizeInds(QAQC_gtype)
write.csv(smry_ind, paste(description, "_", date, "_ind_summary.csv", sep = "",
                          row.names = FALSE))


dups <-dupGenotypes(QAQC_gtype, num.shared = 0.85)
write.csv(dups, paste(description, "_", date, "_dups.csv", sep = " ",
                      row.names = FALSE))
