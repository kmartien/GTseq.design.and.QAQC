library(ggplot2)
library(tidyverse)
library(dplyr)
library(Mnov.GTseq.data)
library(swfscMisc)
source("R/3_fastq.handling/Count.reads.in.fastq.R")
source("R/functions/summarize.depth.at.target.SNPs.R")

project <- "RunMS58"
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
####################################
# Open terminal window, navigate to folder with bam files, then paste:

for FILE in *.bam; do samtools depth -b /Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.data/data-raw/Mnov.targetSNP.bed $FILE > ${FILE%.bam}.coverage; done

####################################
cov.files <- list.files(path = bam.dir, pattern = "*.coverage")
depth.sum <- summarize.depth.at.target.SNP(bam.dir, cov.files, loci)

depth.sum$depth.per.ind$mean <-  nrow(loci) * depth.sum$depth.per.ind$mean
names(depth.sum$depth.per.ind) <- c("Sample","tot.aligned.reads","median.aligned.reads.per.loc","gt.10.aligned", "gt.20.aligned")
depth.sum$depth.per.ind$Sample <- paste0(depth.sum$depth.per.ind$Sample, "_R1.fastq.gz")

ind.summary <- right_join(plate.info, tot.reads.fastq) %>% relocate(Sample) %>%
  left_join(depth.sum$depth.per.ind)
ind.summary <-  left_join(tot.reads.fastq, depth.sum$depth.per.ind)
ind.summary$pct.on.target <-  ind.summary$tot.aligned.reads/ind.summary$tot.reads
loc.summary <- depth.sum$loc.depth

# save results
write.csv(ind.summary, file = paste0("results-raw/", project, ".ind.summary.csv"))
write.csv(loc.summary, file = paste0("results-raw/", project, ".loc.summary.csv"))
save(ind.summary, loc.summary, file = paste0("results-R/", project, ".fastq.BAM.summary.rda"))

# plot results

line.labels <- data.frame(x = c(20,20), y = c(300,200), label = c("80% of loci", "50% of loci"))
g.num.genotypable <- ggplot(ind.summary) + geom_bar(aes(x = reorder(Sample, gt.10.aligned), y = gt.10.aligned), stat = "identity") +
  geom_hline(yintercept = 294) + geom_hline(yintercept = 184) +
  labs(title = project, x = "Ranked samples", y = "Target SNPs with >=10 reads") +
  geom_text(data = line.labels, aes(x = x, y = y, label = label)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

pdf(file = paste0("results-raw/", project, ".summary.pdf"))
depth.sum$plots
g.num.genotypable
dev.off()  
