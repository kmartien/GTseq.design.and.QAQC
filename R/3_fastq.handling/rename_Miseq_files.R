library(swfscMisc)
library(dplyr)

fastq.path <- "/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/Data/Fastq/RunMS43"
fnames <- data.frame(list.files(path = paste0(fastq.path,"/OG files"),
                                pattern=".fastq"))
names(fnames) <- "og.fnames"

fname.parts <- do.call('rbind', lapply(strsplit(fnames$og.fnames, split = "_"), function(i){i}))
LABID <- do.call('c', strsplit(fname.parts[,1], split = "b")) %>% as.numeric()
fname.parts[,1] <- paste0("z0", zero.pad(LABID))

new.names <- sapply(1:dim(fname.parts)[1], function(i){
  paste(fname.parts[i,1], fname.parts[i,2], fname.parts[i,3], fname.parts[i,4], fname.parts[i,5], sep="_")
#  do.call(paste, lapply(1:dim(fname.parts)[2], function(j){fname.parts[i,j]}))
})

fnames$new.names <- new.names
for (i in 1:nrow(fnames)){
  file.rename(from = paste0(fastq.path,"/",fnames$og.fnames[i]), to = paste0(fastq.path, "/", fnames$new.names[i]))
}
