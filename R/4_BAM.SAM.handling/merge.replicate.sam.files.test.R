library(dplyr)

project <- "merge.sam.test"

label.file <- read.table(paste0("data-raw/mplot_labels/", project, ".label.txt"))

counts <- group_by(label.file, V2) %>% summarise(length(V1))

replicates <- counts$V2[which(counts[,2] > 1)]

rep.files <- filter(label.file, V2 %in% replicates)

## Run this code, then take all of the comments it prints out and copy and paste them into Terminal (one at a time)
new.label.file.rows <- do.call(rbind, lapply(1:length(replicates), function(i){
  f <- rep.files$V1[which(rep.files$V2 == replicates[i])]
  out.f <- paste0(replicates[i], "_merged.sam")
  system(paste("samtools merge", out.f, f[1], f[2], sep=" "))
  system2(command = "samtools",
          args = c("merge", paste0("data-raw/sam.files/", project, "/", replicates[i], "_merged.sam"),
                   paste0("data-raw/sam.files/", project, "/", f[1]), 
                   paste0("data-raw/sam.files/", project, "/", f[2])), 
          stdout = TRUE, stderr = TRUE)
  c(out.f, replicates[i], NA)
}))

write.table(label.file, file = paste0("microhaplot/", project, ".replicates.separate.label.txt"),
                                      col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")

label.file <- label.file[-which(label.file$V2 %in% replicates),]
label.file <- rbind(label.file, new.label.file.rows)
write.table(label.file, file = paste0("microhaplot/", project, ".label.txt"), 
                                      col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
