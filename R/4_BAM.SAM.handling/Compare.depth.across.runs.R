library(tidyverse)
library(ggplot2)
library(grid)
library(gridExtra)

seq.projects <- c("GTseq.val.all", "GTseq.prod.all", "RunMS43")
fastq.summaries <- c("HB_val_fastqs", "HBW_Production_fastqs", "RunMS43")


#depth.long <- do.call(rbind, lapply(1:length(seq.projects), function(r){
run.dat <- do.call(rbind, lapply(1:length(seq.projects), function(r){
  load(paste0("results-R/", seq.projects[r], ".depth.at.target.SNPs.rda"))
  dpth <- select(loc.depth, c(locus, mean)) %>% bind_cols(seq.run = seq.projects[r])
  load(paste0("results-R/Fastq.summary.", fastq.summaries[r], ".primer.sequences.rda"))
  loc.sum <- left_join(dpth, sum.by.loc)
}))
names(run.dat) <- c("locus", "mean.depth", "seq.run", "fwd.reads", "rev.reads", "read.ratio")
run.dat$read.ratio[which(run.dat$read.ratio == Inf)] <- NaN

plots <- list()
plots[[1]] <- ggplot(run.dat) + 
  geom_bar(aes(x = reorder(locus, mean.depth), y = mean.depth, fill = seq.run), stat = "identity") +
  theme(legend.position = c(0.25,0.75))
plots[[2]] <- ggplot(filter(run.dat, seq.run %in% seq.projects[1:2])) + 
  geom_bar(aes(x = reorder(locus, mean.depth), y = read.ratio, fill = seq.run), stat = "identity", position = "dodge") +
  theme(legend.position = c(0.25,0.75))
plots[[3]] <- ggplot(filter(run.dat, seq.run %in% seq.projects[1:2])) + 
  geom_bar(aes(x = reorder(locus, mean.depth), y = read.ratio, fill = seq.run), stat = "identity", position = "dodge") +
  theme(legend.position = c(0.25,0.75)) + ylim(c(0,100))

depth.wide <- select(run.dat, c(locus, mean.depth, seq.run)) %>% pivot_wider(names_from = seq.run, values_from = mean.depth)

plots[[4]] <- ggplot(depth.wide, aes(x = GTseq.val.all, y = GTseq.prod.all)) + 
  geom_point() +   labs(title = "Mean Read Depth") 
plots[[5]] <- ggplot(depth.wide, aes(x = GTseq.val.all, y = RunMS43)) + geom_point()
plots[[6]] <- ggplot(depth.wide, aes(x = GTseq.prod.all, y = RunMS43)) + geom_point()

plots$nrow <- 3
pdf("results-raw/run.depth.comparison.pdf")
do.call(grid.arrange, plots)
dev.off()

ratio.wide <- select(run.dat, c(locus, read.ratio, seq.run)) %>% pivot_wider(names_from = seq.run, values_from = read.ratio)
GTSEEK.primer.test <- read.csv("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.gtseq.data/data-raw/GTseq_val-primer_test.csv")
load("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/Mnov.GTseq.data/data/primer.names.rda")
primer.names$LocusID <- paste0(primer.names$CHROM, "_", primer.names$POS)
names(primer.names)[3] <- "locus"

GTSEEK.primer.test <- left_join(GTSEEK.primer.test, select(primer.names, c(locus, LocusID)))
ratio.wide <- left_join(select(ratio.wide, -RunMS43), select(GTSEEK.primer.test, -LocusID))
ratio.wide$GTSEEK.ratio <- ratio.wide$REV / ratio.wide$FWD
for (i in 1:nrow(ratio.wide)) {if(ratio.wide$GTSEEK.ratio[i] < 1) ratio.wide$GTSEEK.ratio[i] <- 1/ratio.wide$GTSEEK.ratio[i]}

save(run.dat, depth.wide, ratio.wide, file = "data/compare.depth.across.runs.rda")
