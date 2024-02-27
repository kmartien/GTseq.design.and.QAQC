# Use bash script bash/merge.sam.files.sh to merge the files. Use this script to
# generate a microhaplot label file that uses the merged files for replicated
# individuals

library(dplyr)
library(tidyverse)

sam.folder <- "RunMS58"

all.files <- data.frame(fname = list.files(path = paste0("data-raw/sam.files/", sam.folder), pattern = ".sam"))
all.files$LABIDs <- substr(all.files$fname, start = 1, stop = 8)
labids <- unique(all.files$LABIDs)

mplot.labels <- do.call('rbind', lapply(labids, function(i){
  file.list <- filter(all.files, LABIDs %in% i)
  if (nrow(file.list) == 1) {return(file.list)} else {
    return(file.list[grep("merged", file.list$fname),])
  }
}))
mplot.labels$stratum <- NA

write.table(mplot.labels, file = paste0("data-raw/mplot_labels/", sam.folder, ".label.txt"), 
                                      col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
