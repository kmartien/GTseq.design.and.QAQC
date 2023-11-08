# NOTE: This script assumes you've already samtools depth command to calculate
# coverage from BAM files, e.g., by executing this line in a terminal window:
# for FILE in *.bam; do samtools depth -b /Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.data/data-raw/Mnov.targetSNP.bed $FILE > ${FILE%.bam}.coverage; done

library(ggplot2)
library(tidyverse)
library(dplyr)
library(Mnov.GTseq.data)
library(swfscMisc)
source("R/3_fastq.handling/Count.reads.in.fastq.R")
source("R/functions/summarize.depth.at.target.SNPs.R")

project <- "RunMS51"
read.type <- "paired" # single or paired
fastq.dir <- paste0("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/Data/Fastq/", project)
bam.dir <- paste0("data-raw/bam.files/", project)
data("SWFSC.bed")

# calculate total reads in each fastq
fastq.files <- list.files(path = fastq.dir, pattern = "*.fastq.gz")
if(read.type == "paired") {
  fastq.files <- fastq.files[grep("R1", fastq.files)]
}
tot.reads.fastq <- count.reads.in.fastq(fastq.dir, fastq.files)

# summarize depth at target SNPs for BAM alignment files
loci <- SWFSC.bed
names(loci) <- c("locus", "start", "stop")
cov.files <- list.files(path = bam.dir, pattern = "*.coverage")
depth.sum <- summarize.depth.at.target.SNP(bam.dir, cov.files, loci)

# combine plate info with per individual summaries
plate.info <- read.csv("data-raw/Mnov2023_Prod_plate-map.csv")
names(plate.info)[6] <- "Sample"
depth.sum$depth.per.ind$mean <-  nrow(loci) * depth.sum$depth.per.ind$mean
names(depth.sum$depth.per.ind) <- c("Sample","tot.aligned.reads","median.aligned.reads.per.loc","gt.10.aligned", "gt.20.aligned")
depth.sum$depth.per.ind$Sample <- paste0(depth.sum$depth.per.ind$Sample, "_R1.fastq.gz")

RunMS51.ind.summary <- right_join(plate.info, tot.reads.fastq) %>% relocate(Sample) %>%
  left_join(depth.sum$depth.per.ind)
RunMS51.ind.summary$pct.on.target <-  RunMS51.ind.summary$tot.aligned.reads/RunMS51.ind.summary$tot.reads
RunMS51.loc.summary <- depth.sum$loc.depth

# save results
write.csv(RunMS51.ind.summary, file = "results-raw/RunMS51.ind.summary.csv")
write.csv(RunMS51.loc.summary, file = "results-raw/RunMS51.loc.summary.csv")
save(RunMS51.ind.summary, RunMS51.loc.summary, file = "results-R/RunMS51.fastq.BAM.summary.rda")

# plot results

RunMS51.ind.summary$well <- paste0(RunMS51.ind.summary$Row, RunMS51.ind.summary$Column)
line.labels <- data.frame(x = c(20,20), y = c(300,200), label = c("80% of loci", "50% of loci"))
g.num.genotypable <- ggplot(RunMS51.ind.summary) + geom_bar(aes(x = reorder(Miseq.num, gt.10.aligned), y = gt.10.aligned), stat = "identity") +
  geom_hline(yintercept = 294) + geom_hline(yintercept = 184) +
  labs(title = "RunMS51", x = "Ranked samples", y = "Target SNPs with >=10 reads") +
  geom_text(data = line.labels, aes(x = x, y = y, label = label)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

plots <- lapply(c(1,2,3,8), function (p){
  dat <- filter(RunMS51.ind.summary, Plate == p)
  ggplot(dat) + geom_bar(aes(x = reorder(well, Miseq.num), y = gt.10.aligned), stat = "identity") +
    geom_hline(yintercept = 294) + geom_hline(yintercept = 184) +
    labs(title = paste0("RunMS51 plate ", p), x = "Loci ordered by column", y = "Target SNPs with >=10 reads") +
    scale_x_discrete(breaks = paste0("A", 1:(length(dat$well)/8)))
})

pdf(file = paste0("results-raw/RunMS51.summary.pdf"))
depth.sum$plots
g.num.genotypable
plots
dev.off()  
