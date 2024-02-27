library(dplyr)
library(Mnov.GTseq.data)
source("R/3_fastq.handling/Count.reads.in.fastq.R")
source("R/3_fastq.handling/Find.primer.seqs.in.fastq.R")

project <- "RunMS58"
read.type <- "paired" # single or paired
fastq.dir <- paste0("/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/Data/Fastq/", project)

#primer.file <- "Illumina.old.smallRNA.primer"
data("SWFSC.bed")
primer.file <- "primer.sequences"
primers <- read.csv(paste0("data-raw/", primer.file, ".csv")) %>% 
  filter(locus %in% SWFSC.bed$locus)

fastq.files <- list.files(path = fastq.dir, pattern = "*.fastq.gz")

# use these lines to limit the number of files and primers for testing purposes
#fastq.files <- fastq.files[97:200]
primers <- primers[1:10,]

if(read.type == "paired") {
  fastq.files <- fastq.files[grep("R1", fastq.files)]
}

date()
all.matches <- primer.seqs.in.fastq(fastq.dir, fastq.files, primers, read.type)
date()

tot.reads <- count.reads.in.fastq(fastq.dir, fastq.files)

rev <- all.matches %>% group_by(Sample) %>% summarise(rev = sum(rev))
sum.by.ind <- all.matches %>% group_by(Sample) %>% summarise(fwd = sum(fwd)) %>%
  left_join(rev) %>% left_join(tot.reads)

rev <- all.matches %>% group_by(locus) %>% summarise(rev = sum(rev))
sum.by.loc <- all.matches %>% group_by(locus) %>% summarise(fwd = sum(fwd)) %>%
  left_join(rev)
sum.by.loc$ratio <- sapply(1:nrow(sum.by.loc), function(l){
  max(sum.by.loc[l,2:3]) / min(sum.by.loc[l,2:3])
})

save(all.matches, sum.by.ind, sum.by.loc, file = paste0("results-R/Fastq.summary.", project, ".", primer.file,".rda"))

