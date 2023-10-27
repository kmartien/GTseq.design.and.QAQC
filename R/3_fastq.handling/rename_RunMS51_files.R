library(swfscMisc)
library(dplyr)

fastq.path <- "/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/Data/Fastq/RunMS51"
fnames <- data.frame(list.files(path = paste0(fastq.path,"/OG files"),
                                pattern=".fastq"))
names(fnames) <- "og.fnames"
plate.info <- read.csv("data-raw/Mnov2023_Prod_plate-map.csv")

fname.parts <- do.call('rbind', lapply(strsplit(fnames$og.fnames, split = "_"), function(i){i}))
Miseq.num <- as.numeric(do.call('rbind', lapply(strsplit(fname.parts[,4], split = "n00"), function(i){i[2]})))
Miseq.num[which(fname.parts[,5] == "R2.fastq.gz")] <- Miseq.num[which(fname.parts[,5] == "R2.fastq.gz")] - 336
fname.parts[,4] <- paste0("n00", as.character(Miseq.num))

fname.parts[,3] <- "MS51"

new.names <- sapply(1:dim(fname.parts)[1], function(i){
  paste(fname.parts[i,1], fname.parts[i,2], fname.parts[i,3], fname.parts[i,4], fname.parts[i,5], sep="_")
})

fnames$new.names <- new.names
for (i in 1:nrow(fnames)){
  file.rename(from = paste0(fastq.path,"/",fnames$og.fnames[i]), to = paste0(fastq.path, "/", fnames$new.names[i]))
}
