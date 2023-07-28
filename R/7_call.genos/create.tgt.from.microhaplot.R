library(tidyverse)
library(dplyr)
source("R/functions/tgt.2.geno.table.R")
source("R/functions/haps.w.ns.R")

project <- "RunMs45.90s.528.allSNPs"

AB.min.het <- 3/7
AB.max.homo <- 2/8
min.AR.het <- 3/10
max.AR.homo <- 2/10
min.read.depth <- 10
num.locs <- 384
min.genos.per.ind <- 307

#locus_annotation <- readRDS(paste0("data/", project, "-microhaplot-locus_annotation.rda"))
genos <- readRDS(paste0("/Users/Shared/KKMDocuments/Documents/Github.Repos/Shiny/microhaplot/", project, ".rds"))[,-1]
genos <- mutate(genos, id.loc = paste0(id, "-", locus))
genos$to.remove <- FALSE

# remove haplotypes that are NA
genos$to.remove[which(is.na(genos$haplo))] <- TRUE
length(which(genos$to.remove))

# remove genotypes that don't meet the minimum read depth criterion
inds <- unique(genos$id)
loci <- unique(genos$locus)
tot.depth <- data.frame(do.call(rbind, lapply(inds, function(i){
  df <- filter(genos, id == i) %>% filter(to.remove == FALSE)
  do.call(rbind, lapply(loci, function(l){
    haps <- filter(df, locus == l)
    data.frame(i, l, sum(haps$depth))
  }))
})))
names(tot.depth) <- c("id", "locus", "tot.depth")

genos.to.drop <- filter(tot.depth, tot.depth < min.read.depth) %>% mutate(id.loc = paste0(id, "-", locus))
genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
length(which(genos$to.remove))

# remove genotypes for individuals whose read depth of its major haplotype at a locus is less than min.read.depth/2
genos.to.drop <- filter(genos, depth < min.read.depth/2) %>% filter(rank == 1) %>% filter(to.remove == FALSE)
genos$to.remove[which(genos$id.loc %in% genos.to.drop$id.loc)] <- TRUE
length(which(genos$to.remove))

# remove haplotypes with rank > 1 and AB < AB.min.het 
genos$to.remove[which(genos$rank > 1 & genos$allele.balance < AB.min.het)] <- TRUE
length(which(genos$to.remove))

# check loci for which some individuals have more than two haplotypes, decide which to drop
excess.haps <- filter(genos, rank > 2) %>% filter(to.remove == FALSE)
inds.to.check <- filter(genos, id.loc %in% excess.haps$id.loc)
locs.to.check <- data.frame(table(inds.to.check$locus)) %>% mutate(Status = "Accept")
names(locs.to.check)[1] <- "locus"
#write.csv(locs.to.check, file = paste0("data-raw/", project, ".locs.to.check.csv"), row.names = FALSE)

### SEE IF I CAN FIGURE OUT WHAT'S GOING ON 155, 219, AND 302
#locs.to.drop <- read.csv(file = paste0("data-raw/", project, ".locs.to.check.csv")) %>%
locs.to.drop <- read.csv(file = paste0("data-raw/all.locs.to.check.csv")) %>%
  filter(Status %in% c("Reject", "??")) %>% select(locus)

genos$to.remove[which(genos$locus %in% locs.to.drop$locus)] <- TRUE

# remove genotypes with rank > 4
genos$to.remove[which(genos$rank > 4)] <- TRUE
length(which(genos$to.remove))

# remove filtered genotypes
genos.filtered <- filter(genos, to.remove == FALSE) 
loci <- unique(genos.filtered$locus)

# check loci for which some individuals have more than two haplotypes, decide which to drop
excess.haps <- filter(genos.filtered, rank > 2) %>% filter(to.remove == FALSE)
inds.to.check <- filter(genos.filtered, id.loc %in% excess.haps$id.loc)
locs.to.check <- data.frame(table(inds.to.check$locus)) %>% mutate(Status = "Accept")
names(locs.to.check)[1] <- "locus"

# calculate haplotype freqs at remaining loci
#hapfreqs <- lapply(loci, function(l){
#  df <- filter(genos.filtered, locus == l)
#  hapfreq <- table(df$haplo)
#  return(sort(hapfreq, decreasing = TRUE))
#})
#names(hapfreqs) <- loci

# create tgt-like data structure
tgt <- data.frame(do.call(rbind, lapply(unique(genos.filtered$id.loc), function(x){
  df <- filter(genos.filtered, id.loc == x)
  loc <- df$locus[1]
  ind <- df$id[1]
  haps <- df$haplo
  if(length(haps) < 4) haps <- c(haps, rep(NA, (4-length(haps))))
  depth <- df$depth
  if(length(depth) < 4) depth <- c(depth, rep(0, (4-length(depth))))
  res <- c(loc, ind, haps, depth)
  names(res) <- c("locus", "Indiv", "haplo.1", "haplo.2", "haplo.3", "haplo.4", 
                  "depth.1", "depth.2", "depth.3", "depth.4")
  if(nrow(df) == 1){
    res[4] <- res[3] 
    res[8] <- 0
  }
  return(res)
})))

tgt$depth.1 <- as.integer(tgt$depth.1)
tgt$depth.2 <- as.integer(tgt$depth.2)
tgt$depth.3 <- as.integer(tgt$depth.3)
tgt$depth.4 <- as.integer(tgt$depth.4)
#tgt$allele.balance <- tgt$depth.1/tgt$depth.2

# summarize missing data
missing.data.ind <- data.frame(table(tgt$Indiv)) %>%
  mutate(missing = num.locs-Freq)
names(missing.data.ind) <- c("labID", "genos", "missing")
length(which(missing.data.ind$genos >= min.genos.per.ind))
rejected.inds <- missing.data.ind$labID[which(missing.data.ind$genos < min.genos.per.ind)]

tgt <- tgt[-which(tgt$Indiv %in% rejected.inds),]
num.inds <- length(unique(tgt$Indiv))

missing.data.loc <- data.frame(table(tgt$locus)) %>%
  mutate(missing = num.inds-Freq)
names(missing.data.loc)[1] <- "locus"
length(which(missing.data.loc$missing > round(num.inds * 0.2)))

geno.table <- tgt.2.geno.table(tgt)

save(geno.table, tgt, genos.filtered, file = paste0("results-R/", project, ".", min.read.depth, "readsMin.tgt.rda"))
write.csv(missing.data.ind, file = paste0("results-raw/", project, ".", min.read.depth, "readsMin.num.genos.per.ind.csv"))

