library(ggplot2)
library(grid)
library(gridExtra)

#source("R/SNPFILT/libraries.R")
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/SNPFILT/ggplot.R")
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/SNPFILT/VCFfilterstats.R")
source("/Users/Shared/KKMDocuments/Documents/Github.Repos/Mnov/GTseq.design.and.QAQC/R/SNPFILT/xtrafunctions.R")

visualize.vcf.ind <- function(res.dir, fname){

  pdfname = paste(res.dir,"/",fname,"_summary_Indivs.pdf",sep="")
  pdf(pdfname) 
  
  ind_stats_raw <- read.ind.stats(dir = res.dir, fname)
  loc_stats_raw <- read.loc.stats(dir = res.dir, fname) 
  names(ind_stats_raw) <- sapply(names(ind_stats_raw), function(x) strsplit(x,fname))
  names(loc_stats_raw) <- sapply(names(loc_stats_raw), function(x) strsplit(x,fname))
  
  # plot missing data per indv ----
  p1 <- ggplot(ind_stats_raw, aes(x = MISS_)) +
    geom_histogram(binwidth = .01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MISS_, 
                                     na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 0.5), color = "darkblue", 
               linetype = "dashed",  linewidth = 1) + 
    labs(title = paste("# individual = ", dim(ind_stats_raw)[1],sep=""), x = "missing data per indv") #+ theme_standard
  
  # plot Fis per indv ----
  p2 <- ggplot(ind_stats_raw, aes(x = Fis_)) +
    geom_histogram(binwidth = .01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(Fis_, 
                                     na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 0), color = "darkblue", linetype = "dashed", 
               linewidth = 1) + labs(x = "Fis per indv") #+
    #theme_standard
  
  # plot read depth per indv ----
  p3 <- ggplot(ind_stats_raw, aes(x = MEAN_DEPTH_)) +
    geom_histogram(binwidth = 10, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH_, 
                                     na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 20),
               color = "darkblue", linetype = "dashed", linewidth = 1) +
    labs(x = "mean read depth per indv") #+ theme_standard
  
  # plot depth vs missing ----
  p4 <- ggplot(ind_stats_raw, aes(x = MEAN_DEPTH_, y = MISS_)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH_, 
                                     na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 20), color = "darkblue", 
               linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = mean(MISS_, 
                                     na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = 0.5), color = "darkblue", 
               linetype = "dashed", linewidth = 1) +
    labs(x = "mean depth per indv", y = "% missing data") #+
    #theme_standard
  
  # plot Fis vs missing data per indv ----
  p5 <- ggplot(ind_stats_raw, aes(x = Fis_, 
                                  y =  MISS_)) + geom_point() +
    geom_vline(aes(xintercept = mean(Fis_, 
                                     na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 0), color = "darkblue", linetype = "dashed", 
               linewidth = 1) + geom_hline(aes(yintercept = mean(MISS_,   na.rm = TRUE)),
                                      color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = 0.5),
               color = "darkblue", linetype = "dashed", linewidth = 1) +
    labs(x = "Fis per indv", y = "% missing data") #+
#    theme_standard
  
  # plot Fis vs mean depth per indv ----
  p6 <- ggplot(ind_stats_raw, aes(x = Fis_, 
                                  y =  MEAN_DEPTH_)) +
    geom_point() +
    geom_vline(aes(xintercept = mean(Fis_, 
                                     na.rm = TRUE)),
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 0),
               color = "darkblue", linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = mean(MEAN_DEPTH_, 
                                     na.rm = TRUE)),
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = 20),
               color = "darkblue", linetype = "dashed", linewidth = 1) +
    labs(x = "Fis per indv", y = "mean depth per indv") #+
#    theme_standard
  
#  ggdraw() + draw_plot(p1) + draw_plot(p2) + draw_plot(p3) + draw_plot(p4) + 
#    draw_plot(p5) + draw_plot(p6)
#  m1 <- multiplot(p1, p2, p3, p4, p5, p6, cols=2)
  plots <- list(p1, p2, p3, p4, p5, p6)
  plots$nrow <- 3
  do.call(grid.arrange, plots)
  dev.off()
}

# VISUALIZE LOCUS SUMMARY STATS
visualize.vcf.loci <- function(res.dir, fname){
  
  pdfname = paste(res.dir,"/",fname,"_summary_Loci.pdf",sep="")
  pdf(pdfname) 
  
  ind_stats_raw <- read.ind.stats(dir = res.dir, fname)
  loc_stats_raw <- read.loc.stats(dir = res.dir, fname) 
  names(ind_stats_raw) <- sapply(names(ind_stats_raw), function(x) strsplit(x,fname))
  names(loc_stats_raw) <- sapply(names(loc_stats_raw), function(x) strsplit(x,fname))
  
  # plot distribution missing data per locus ----
  p7 <- ggplot(loc_stats_raw, aes(x = MISS_)) +
    geom_histogram(binwidth = 0.01, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MISS_, na.rm = TRUE)),color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 0.1), color = "darkblue", linetype = "dashed", linewidth = 1) +
    labs(title = paste("# loci = ", dim(loc_stats_raw)[1],sep=""), x = "% missing data per locus") #+
#    theme_standard
  
  # plot distribution mean read depth ----
  lmiss.95pctile <- sort(loc_stats_raw$MEAN_DEPTH_)[.95*length(loc_stats_raw$MEAN_DEPTH_)]
  p8 <- ggplot(loc_stats_raw, aes(x = MEAN_DEPTH_)) +
    geom_histogram(binwidth = 5, color = "black", fill = "grey95") +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH_, na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 20), color = "darkblue", linetype = "dashed", linewidth = 1) + 
    labs(title = paste("95th percentile meanDP = ",lmiss.95pctile,sep=""), x = "mean read depth per locus") #+
#    theme_standard
  
  # plot read depth vs missing data ----
  p9 <- ggplot(loc_stats_raw, aes(x = MEAN_DEPTH_, y = MISS_)) + geom_point() +
    geom_vline(aes(xintercept = mean(MEAN_DEPTH_, na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = 20), color = "darkblue", linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = mean(MISS_, na.rm = TRUE)), color = "red", linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = 0.1), color = "darkblue", linetype = "dashed", linewidth = 1) +
    labs(x = "mean depth per locus", y = "% missing data") #+
#    theme_standard()
  
  # plot no of SNPs per contig ----
  p10 <- loc_stats_raw %>%
    dplyr::group_by(CHR) %>% dplyr::summarise(Count = dplyr::n()) %>%
    #count("CHR") %>%
    ggplot(aes(x = Count)) +
    geom_histogram(binwidth = 1, color = "black", fill = "grey95") + 
    labs(x = "number of SNPs per contig") #+
#    theme_standard
  temp <- loc_stats_raw %>%
    #count("CHR")
    dplyr::group_by(CHR) %>% dplyr::summarise(Count = dplyr::n())
  
  # plot number of SNPs per contig vs. mean depth ----
  p11 <- left_join(temp, loc_stats_raw, multiple = "all") %>%
    ggplot() +
    geom_point(aes(x = Count, y = MEAN_DEPTH_)) +
    labs(x = "number of SNPs per contig", y = "mean depth") #+
#    theme_standard
  
  # plot depth vs SNP quality ----
  site_qual <- read.table(paste(results.dir, "/", fname,".lqual",sep=""), 
                          header = TRUE, stringsAsFactors = FALSE)  %>%
    mutate(PROB = 10^(-QUAL/10))
  temp <- data.frame(loc_stats_raw$MEAN_DEPTH_, site_qual$QUAL)
  colnames(temp)<-c("depth", "qual")
  high_depth <- mean(loc_stats_raw$MEAN_DEPTH) + sqrt(var(loc_stats_raw$MEAN_DEPTH))

  p12 <- ggplot(temp, aes(x = depth, y = qual)) +
    geom_point(size = 1) +
    geom_vline(aes(xintercept = mean(depth, na.rm = TRUE)),
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(aes(xintercept = high_depth),
               color = "darkblue", linetype = "dashed", linewidth = 1) +
    geom_hline(aes(yintercept = mean(qual, na.rm = TRUE)),
               color = "red", linetype = "dashed", linewidth = 1) +
    geom_abline(aes(slope = 2,intercept = 0),
               color = "darkblue", linetype = "dashed", linewidth = 1) +
    labs(x = "mean depth per locus", y = "SNP quality") #+
#    theme_standard
#  m2 <- multiplot(p7, p8, p9, p10, p11, p12, cols=2)
  plots <- list(p7, p8, p9, p10, p11, p12)
  plots$nrow <- 3
  do.call(grid.arrange, plots)
  dev.off()
}