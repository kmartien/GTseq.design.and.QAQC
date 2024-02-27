library(dplyr)
library(ggplot2)
library(gridExtra)
library(Mnov.GTseq.data)
library(swfscMisc)

data("GTseq.samps.final")
data("tissue.archive.data")

projects <- list("GTseq.val","RunMS51","RunMS54")

CCES.samps <- filter(GTseq.samps.final, USE == "CCES") %>% select(c(LABID, Duplicate1))
CCES.samps <- c(CCES.samps$LABID, CCES.samps$Duplicate1) %>% na.omit()
CCES.samps <- paste0("z0", zero.pad(as.numeric(CCES.samps)))
tissue.archive.data$LABID <- paste0("z0", zero.pad(tissue.archive.data$LABID))

diluted.samples <- read.csv("data-raw/diluted.samples.csv") %>% select(c(LAB.ID,dil))
names(diluted.samples) <- c("LABID","ratio")
diluted.samples$LABID <- paste0("z0", zero.pad((diluted.samples$LABID)))

tgt.list <- lapply(projects, function(p){
  load(paste0("results-R/", p, ".20readsMin.geno.eval.rda"))
  tgt$Run <- p
  return(tgt)
})
names(tgt.list) <- projects

tgt <- do.call(rbind, tgt.list) %>% filter(Indiv %in% CCES.samps)

CCES.dat <- do.call(rbind, lapply(projects, function(p){
  x <- filter(tgt, Run == p)
  missing.data.ind <- data.frame(table(x$Indiv))
  names(missing.data.ind) <- c("LABID","genos")
  dat <- left_join(missing.data.ind, tissue.archive.data) %>% select(c(LABID, genos, YR, MO, DA))
  dat$date <- as.Date(paste(dat$YR, zero.pad(dat$MO), zero.pad(dat$DA), sep = "-"))
  dat$julian.date <- format(dat$date, "%j")
  dat$Run <- p
  return(dat)
})) %>% arrange(LABID)
CCES.dat$diluted <- "No"
CCES.dat$diluted[which(CCES.dat$LABID %in% diluted.samples$LABID)] <- "Yes"
CCES.dat$diluted[which(CCES.dat$Run == "GTseq.val")] <- "Yes.GTSEEK"
CCES.dat$size <- 3
CCES.dat$size[which(CCES.dat$Run == "Gtseq.val")] <- 6


x.labs <- c("2018-07-04", "2018-07-17", "2018-08-02", "2018-08-10", "2018-09-08", 
            "2018-09-21","2018-10-22")

g.date <-   ggplot(data = CCES.dat, aes(x = julian.date, y = genos, colour = Run)) + 
  geom_point(size = 3, position = "jitter") + scale_x_discrete(labels = x.labs, breaks = format(as.Date(x.labs), "%j")) +
  labs(title = "All CCES", x = "Collection date", y = "Num. loci genotyped") +
  theme(text = element_text(size = 20))  
tiff("results-raw/CCES.by.date.tiff", height = 900, width = 1800)
g.date
dev.off()

x.labs <- CCES.dat$LABID[seq(1, 339, by = 33)]
g.labid <-   ggplot(data = CCES.dat, aes(x = LABID, y = genos, colour = diluted)) + 
  geom_point(size = CCES.dat$size, position = "jitter") + scale_x_discrete(labels = x.labs, breaks = x.labs) +
  geom_point(data = filter(CCES.dat, Run == "GTseq.val"), aes(x = LABID, y = genos), size = 6, colour = "blue") +
  geom_point(data = filter(CCES.dat, diluted == "Yes"), aes(x = LABID, y = Concentration), size = 6, colour = "green") +
  labs(title = "All CCES", x = "LABID", y = "Num. loci genotyped") +
  theme(text = element_text(size = 30))  

Miseq.samp.summary <- read.csv("results-raw/MiSeq.samples.summary.csv")
CCES.dat <- left_join(CCES.dat, select(Miseq.samp.summary, c(LABID,Concentration)))

g.concentration <-   ggplot(data = CCES.dat, aes(x = LABID, y = Concentration, colour = diluted)) + 
  geom_point(size = 3, position = "jitter") + scale_x_discrete(labels = x.labs, breaks = x.labs) +
  geom_point(data = filter(CCES.dat, Run == "GTseq.val"), aes(x = LABID, y = Concentration), size = 6, colour = "blue") +
  geom_point(data = filter(CCES.dat, diluted == "Yes"), aes(x = LABID, y = Concentration), size = 6, colour = "green") +
  labs(title = "All CCES", x = "LABID", y = "Concentration") +
  theme(text = element_text(size = 30))  

plots <- list(g.labid, g.concentration)
plots$nrow <- 2
tiff("results-raw/CCES.by.labid.dilution.tiff", height = 1800, width = 1800)
do.call(grid.arrange, plots)
dev.off()

mean.num.genos <- group_by(CCES.dat, diluted) %>% summarise(mean = mean(genos))

