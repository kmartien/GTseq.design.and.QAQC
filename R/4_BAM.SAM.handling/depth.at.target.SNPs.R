# This script uses samtools to calculated read depth for each individual at each locus,
# then calculates summary statistics and produces summary graphs.
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

#data("SWFSC.bed")
project <- "RunMS45.528locs"
which.files <- "90s" # "all", "30s", or "90s"
bed.file <- "/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.data/data-raw/528locs.allSNPs.bed"

####################################
# Open terminal window, navigate to folder with bam files, then paste:

for FILE in *.bam; do samtools depth -b /Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.data/data-raw/Mnov.targetSNP.bed $FILE > ${FILE%.bam}.coverage; done

####################################

loci <- read.table(bed.file)
#loci <- SWFSC.bed
names(loci) <- c("locus", "start", "stop")

cov.files <- list.files(path = paste0("data-raw/bam.files/", project), pattern = "*.coverage")
if (which.files == "30s") cov.files <- cov.files[seq(from = 1, to = length(cov.files), by = 2)]
if (which.files == "90s") cov.files <- cov.files[seq(from = 2, to = length(cov.files), by = 2)]

cov.list <- lapply(1:length(cov.files), function(i){
  print(i)
#  sample.name <- strsplit(cov.files[i], split = "_")[[1]][1]
  res <- read.table(paste0("data-raw/bam.files/", project, "/", cov.files[i]))
#  names(res) <- c("locus", "POS", paste0("depth.", sample.name))
  names(res) <- c("locus", "POS", cov.files[i])
  return(res[,c(1,3)])
})

depth <- select(loci, c(locus, stop))
for (i in 1:length(cov.list)) {
  depth <- left_join(depth, cov.list[[i]], by = "locus")
}

# replace NAs with zeroes
depth[is.na(depth)] <- 0
depth <- depth[,-2]

# calc mean and median per locus
depth <- rowwise(depth) %>%
  mutate(mean = mean((c_across(where(is.numeric))))) %>%
  mutate(med = median((c_across(where(is.numeric))))) %>%
  arrange(mean)

line.labels <- data.frame(x = c(20,20), y = c(15,25), label = c("10 reads", "20 reads"))

g.mean <- ggplot(depth) + geom_bar(aes(x = reorder(locus, mean), y = mean), stat = "identity") +
  geom_hline(yintercept = 20) + geom_hline(yintercept = 10) +
  labs(title = project, x = "Ranked loci", y = "Mean Coverage at Target SNP") +
  geom_text(data = line.labels, aes(x = x, y = y, label = label)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
g.med <- ggplot(depth) + geom_bar(aes(x = reorder(locus, med), y = med), stat = "identity") +
  geom_hline(yintercept = 20) + geom_hline(yintercept = 10) +
  labs(title = project, x = "Ranked loci", y = "Median Coverage at Target SNP") +
  geom_text(data = line.labels, aes(x = x, y = y, label = label)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Summarize by individual
depth.per.ind <- data.frame(colMeans(depth[,2:(ncol(depth)-2)])) %>%
  bind_cols(apply(depth[,2:(ncol(depth)-2)], 2, median)) %>% rownames_to_column() 
depth.per.ind$rowname <- do.call(rbind, lapply(1:nrow(depth.per.ind), function(i){
  strsplit(depth.per.ind$rowname[i], split = ".f")[[1]][1]
}))

depth.per.ind$gt.10 <- do.call(rbind, lapply(2:(ncol(depth)-2), function(i){
  length(which(depth[,i] > 9))
}))
depth.per.ind$gt.20 <- do.call(rbind, lapply(2:(ncol(depth)-2), function(i){
  length(which(depth[,i] > 19))
}))

names(depth.per.ind) <- c("new.fname", "mean", "median","gt.10", "gt.20")

### THIS BIT IS SPECIFIC TO KAREN'S GTSEEK DATA SUMMARIES
#if (project %in% c("GTseq.val.all", "GTseq.prod.all")) {
#data("gtseq.prod.summary")
#data("gtseq.val.genos")

#gtseq.val.genos$Sample <- do.call(rbind, lapply(1:nrow(gtseq.val.genos), function(i){
#  strsplit(gtseq.val.genos$Sample[i], split = "_p")[[1]][1]
#}))

#load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.GTseq.sandbox/data/labid.to.filename.rda")
#gtseq.smry.all <- rbind(gtseq.prod.summary, gtseq.val.genos[,1:7]) %>%
#  left_join(select(labid.to.filename, c("Sample", "new.fname")))
#gtseq.smry.all$LABID <- paste0("z0", zero.pad(gtseq.smry.all$LABID))
#labid.to.filename$new.fname <- do.call(rbind, lapply(1:nrow(labid.to.filename), function(i){
#  strsplit(labid.to.filename$new.fname[i], split = "_R")[[1]][1]
#}))
#labid.to.filename$Sample <- do.call(rbind, lapply(1:nrow(labid.to.filename), function(i){
#  strsplit(labid.to.filename$Sample[i], split = "_R")[[1]][1]
#}))
#labid.to.filename <- distinct(labid.to.filename)
#names(labid.to.filename) <- c("LABID","new.fname", "Sample", "well","plate")

#depth.per.ind <- left_join(depth.per.ind, select(labid.to.filename, c(new.fname,Sample))) %>%
#  left_join(select(gtseq.smry.all, c(Sample, On.Target.Reads)))
#}
#########################################

num.inds.genotypable.rd10 <- length(which(depth.per.ind$gt.10 > (length(unique(loci$locus)) * 0.8)))
num.inds.genotypable.rd20 <- length(which(depth.per.ind$gt.20 > (length(unique(loci$locus)) * 0.8)))

# SPECIFIC TO KAREN'S GTSEEK SUMMARIES
#g.ind.genos <- ggplot(depth.per.ind, aes(x = On.Target.Reads, y = gt.10, col = "blue")) + 
#  geom_point() + geom_point(aes(x = On.Target.Reads, y = gt.20, col = "red")) +
#  geom_hline(yintercept = 307) + scale_color_manual(labels = c("10", "20"), values = c("blue", "red")) +
#  labs(title = paste(project, "\nNum individuals with read depth >= 10 at 307 loci:", num.inds.genotypable.rd10,
#                     "\nNum individuals with read depth >= 20 at 307 loci:", num.inds.genotypable.rd20),
#                     x = "On Target Reads", y = "Loci Meeting Read Depth Criterion")
g.genos.v.mean.depth <-  ggplot(depth.per.ind, aes(x = mean, y = gt.10, col = "blue")) + 
  geom_point() + geom_point(aes(x = mean, y = gt.20, col = "red")) +
  geom_hline(yintercept = 307) + scale_color_manual(labels = c("10", "20"), values = c("blue", "red")) +
  labs(title = paste(project, "\nNum individuals with read depth >= 10 at 307 loci:", num.inds.genotypable.rd10,
                     "\nNum individuals with read depth >= 20 at 307 loci:", num.inds.genotypable.rd20),
       x = "Mean Read Depth", y = "Loci Meeting Read Depth Criterion")

pdf(file = paste0("results-raw/coverage.summary.",project,".pdf"))
g.mean
g.med
g.genos.v.mean.depth
dev.off()  

loc.depth <- depth
write.csv(loc.depth, file = paste0("results-raw/", project,".",which.files, ".LOCdepth.at.target.SNPs.csv"))
write.csv(depth.per.ind, file = paste0("results-raw/", project,".",which.files, ".INDdepth.at.target.SNPs.csv"))
save(loc.depth, depth.per.ind, file = paste0("results-R/", project,".",which.files, ".depth.at.target.SNPs.rda"))
