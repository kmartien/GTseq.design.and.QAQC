library(ggplot2)
library(dplyr)
library(gridExtra)
library(swfscMisc)

diluted.samples <- read.csv("data-raw/diluted.samples.csv") %>% select(c(LAB.ID,dil))
names(diluted.samples) <- c("LABID","ratio")
diluted.samples$LABID <- paste0("z0", zero.pad((diluted.samples$LABID)))

Miseq.samp.summary <- read.csv("results-raw/MiSeq.samples.summary.csv") %>%
  left_join(diluted.samples)

# all samples, colored by Status
g <- ggplot(Miseq.samp.summary, aes(x = Concentration, y = genos.20readsMin, colour = ratio)) +
  geom_point() +
  scale_x_continuous(breaks=seq(0,1000,100)) +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.9, 0.1)
  )
g

# all samples
g.zoomed <- ggplot(Miseq.samp.summary, aes(x = Concentration, y = genos.20readsMin, colour = Status)) +
  geom_point() + xlim(0,250) +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.9, 0.1)
  )
g.zoomed

g.hist <- ggplot(filter(Miseq.samp.summary, !is.na(Concentration))) + 
  facet_wrap(~Status) +
  geom_histogram(aes(x = Concentration, fill = Status)) +
  xlim(0,150)
g.hist

plots <- list(g, g.zoomed)
plots$nrow <- 2
jpeg("results-raw/Concentration.vs.genotypes.Status.jpg", height = 1000, width = 700)
do.call(grid.arrange, plots)
dev.off()

# all samples, colored by ratio
g <- ggplot(Miseq.samp.summary, aes(x = Concentration, y = genos.20readsMin, colour = ratio)) +
  geom_point() +
  scale_x_continuous(breaks=seq(0,1000,100)) +
  theme(
    text = element_text(size = 25),
    legend.position = c(0.9, 0.1)
  )
g

# zoomed
g.zoomed <- ggplot(Miseq.samp.summary, aes(x = Concentration, y = genos.20readsMin, colour = ratio)) +
  geom_point() + xlim(0,250) +
  theme(
    text = element_text(size = 25),
    legend.position = c(0.9, 0.1)
  )
g.zoomed

plots <- list(g, g.zoomed)
plots$nrow <- 2
jpeg("results-raw/Concentration.vs.genotypes.dilution.ratio.jpg", height = 1000, width = 700)
do.call(grid.arrange, plots)
dev.off()
