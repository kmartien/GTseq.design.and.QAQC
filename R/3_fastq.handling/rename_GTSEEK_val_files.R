library(swfscMisc)
library(dplyr)

fastq.path <- "/Users/Shared/KKMDocuments/Documents/Karen/Structure/Humpbacks/Data/Fastq/HB_val"
fnames <- data.frame(list.files(path = paste0(fastq.path,"/OG files"),
                     pattern=".fastq"))
names(fnames) <- "og.fnames"

load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.Gtseq.data/data/gtseq.val.genos.rda")
gtseq.val.genos <- gtseq.val.genos[-97,1:2]
gtseq.val.genos$LABID <- paste0("z0", zero.pad(gtseq.val.genos$LABID))

base.name.parts <- do.call('rbind', lapply(strsplit(gtseq.val.genos$Sample, split = "_"), function(i){i}))
gtseq.val.genos$base.sample.name <- paste(base.name.parts[,1], base.name.parts[,2], base.name.parts[,3],
                          base.name.parts[,4], base.name.parts[,5], base.name.parts[,6], sep = "_")

fname.parts <- do.call('rbind', lapply(strsplit(fnames$og.fnames, split = "_"), function(i){i}))
fnames$base <- paste(fname.parts[,1], fname.parts[,2], fname.parts[,3],
                         fname.parts[,4], fname.parts[,5], fname.parts[,6], sep="_")

temp <- do.call(rbind, lapply(1:dim(fnames)[1], function(i){
  c(gtseq.val.genos$LABID[which(gtseq.val.genos$base.sample.name == fnames[i,2])],
    paste(gtseq.val.genos$LABID[which(gtseq.val.genos$base.sample.name == fnames[i,2])], fnames[i,1], sep="_"),
    fname.parts[i,6])
}))
fnames <- cbind(fnames, temp)
names(fnames) <- c("Sample", "base", "LABID", "new.fname", "well")
fnames$plate <- "1"

load("data/labid.to.filename.rda")
labid.to.filename$LABID <- paste0("z0", zero.pad(as.numeric(labid.to.filename$LABID)))
labid.to.filename.all <- bind_rows(labid.to.filename, select(fnames, c("LABID", "new.fname", "Sample", "well", "plate")))
labid.to.filename <- labid.to.filename.all
write.csv(labid.to.filename, file = "results/labid.to.filename.csv")
save(labid.to.filename, file = "data/labid.to.filename.rda")
