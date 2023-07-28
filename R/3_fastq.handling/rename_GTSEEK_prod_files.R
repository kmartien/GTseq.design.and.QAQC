library(tidyverse)
library(swfscMisc)

fastq.path <- "/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/Data/Fastq/HBW_Production-fastqs"
fnames <- list.files(path = paste0(fastq.path,"/OG files"),
                     pattern=".fastq")
num.files <- length(fnames)

fname.parts <- do.call('rbind', lapply(strsplit(fnames, split = "_"), function(i){i}))
#labids <- do.call('c', lapply(strsplit(fnames, split = "_"), function(i){i[7]}))
numeric.labids <- as.numeric(fname.parts[,7])
labids <- paste0("z0", zero.pad(numeric.labids))
labids[which(labids == "z0NA")] <- "empty"

new.fnames <- paste(labids, "Mnov_GTseek", fname.parts[,2], fname.parts[,4], "PLT", fname.parts[,6], fname.parts[,8], sep="_")
#file.rename(paste(fastq.path,fnames,sep="/"),paste(fastq.path,new.fnames,sep="/"))

labid.to.filename <- data.frame(cbind(numeric.labids, new.fnames, fnames, fname.parts[,4], fname.parts[,6]))
names(labid.to.filename) <- c("LABID", "new.fname", "Sample", "well", "plate")
write.csv(labid.to.filename, file = "results/labid.to.filename.csv")
save(labid.to.filename, file = "data/labid.to.filename.rda")

plate.info <-distinct(data.frame(cbind(paste(fname.parts[,1], fname.parts[,2], 
                                             fname.parts[,3], fname.parts[,4], 
                                             fname.parts[,5], fname.parts[,6], 
                                             fname.parts[,7],sep= "_"), fname.parts[,4], fname.parts[,6])))
names(plate.info) <- c("Sample", "well", "plate")
write.csv(plate.info, file = "data-raw/plate.info.csv")
save(plate.info, file = "data/plate.info.rds")
