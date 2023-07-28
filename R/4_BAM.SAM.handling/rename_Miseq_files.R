library(swfscMisc)
library(dplyr)

run <- "RunMS45.528locs"
fnames <- data.frame(list.files(path = paste0("data-raw/sam.files/", run, "/"),
                                pattern=".sam"))
names(fnames) <- "og.fnames"

fname.parts <- do.call('rbind', lapply(strsplit(fnames$og.fnames, split = "_"), function(i){i}))
LABID <- do.call('c', strsplit(fname.parts[,1], split = "b")) %>% as.numeric()
fname.parts[,1] <- paste0("z0", zero.pad(LABID))

new.names <- sapply(1:dim(fname.parts)[1], function(i){
  paste(fname.parts[i,1], fname.parts[i,2], fname.parts[i,3], fname.parts[i,4], sep="_")
})

fnames$new.names <- new.names
for (i in 1:nrow(fnames)){
  file.rename(from = paste0("data-raw/sam.files/", run, "/",fnames$og.fnames[i]), to = paste0("data-raw/sam.files/", run, "/", fnames$new.names[i]))
}
