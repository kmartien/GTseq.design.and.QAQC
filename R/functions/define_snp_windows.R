# takes a vcf file and identifies windows of size <= window.sz (default 134bp) 
# containing between min.SNPs and max.SNPs (default 2 -5 SNPs). For each window,
# the position of the first ('start') and last ('stop') SNP in the window is
# specified, along with the total number of SNPs in the window and the mean
# theta_h of those SNPs. Results are returned as a list, with one element per
# chromosome.

# Note that this will crash if your vcf files includes SNPs with more than 
# two alleles. 

define_snp_windows <- function(vcf, window.sz = 134, min.SNPs = 2, max.SNPs = 5){

  # convert vcf to a tibble and add a column 'locus' that combines CHROM and POS
  # of the SNPs
  tidy.vcf <- vcfR2tidy(vcf, single_frame = TRUE, 
                        info_fields = c("DP"), format_fields = c("GT"))$dat %>%
    mutate(locus = paste(CHROM, POS, sep= "_")) %>%
    relocate(locus, .after = POS)

  ## use this line to limit length for testing purposes
  #tidy.vcf <- tidy.vcf[1:100000,]

  ## calculate allele frequencies at each SNP
  allele.freqs <- vcf2allelefreqs(tidy.vcf)

  ## identify unique loci (vcf has entries for every individual at every locus)
  bd <- select(tidy.vcf, c(locus, CHROM, POS)) %>% distinct()
  names(bd) <- c("locus", "CHROM", "start")

  ## identify windows separately on each chromosome
  window.list <- lapply(unique(bd$CHROM), function(l){
    print(l) #printing each locus gives an indication of progress
    bd <- filter(bd, CHROM == l) %>% arrange(start) # make sure SNPs are in order by position

    if (nrow(bd) == 1) return(NULL) #if there's only one SNP on a chromosome, there are no windows
    snp.windows <- data.frame(do.call(bind_rows,lapply(1:(nrow(bd) - 1), function(i){
      ## identify SNPs within window.sz downstream of the focal SNP
      snps.in.window <- filter(bd[(i):(i+5),], start <= bd$start[i] + window.sz)
      
      ## if there are too few or too many SNPs, don't identify a window
      if(nrow(snps.in.window) < min.SNPs | nrow(snps.in.window) > max.SNPs) return(NULL) else {
        
        ## get allele frequencies and calculate theta_h for SNPs in window
        window.af <- filter(allele.freqs, locus %in% snps.in.window$locus)
        h <- mean(sapply(1:nrow(window.af), function(x) {
          theta.h(as.numeric(window.af[x,2:3]))
        }))
        return(bind_cols(bd[i,], stop = max(snps.in.window$start), n.snps = nrow(snps.in.window), mean.theta_h = h))
      }
    })))
    return(snp.windows)
  })
  names(window.list) <- unique(bd$CHROM)
  
  ## remove chromosomes with no windows identified
  window.list <- window.list[-which(sapply(window.list, is.null))]
  return(window.list)
}
