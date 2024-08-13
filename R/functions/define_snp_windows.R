# takes a vcf file and identifies windows of size <= window.sz (default 134bp) 
# containing between min.SNPs and max.SNPs (default 2 -5 SNPs). For each window,
# the position of the first ('start') and last ('stop') SNP in the window is
# specified, along with the total number of SNPs in the window and the mean
# theta_h of those SNPs. Results are returned as a list, with one element per
# chromosome.

# Note that this will crash if your vcf files includes SNPs with more than 
# two alleles. 

define_snp_windows <- function(tidy.vcf, window.sz = 134, min.SNPs = 1, max.SNPs = 5){

  ## calculate allele frequencies at each SNP
  allele.freqs <- vcf2allelefreqs(tidy.vcf)

  ## identify unique loci (vcf has entries for every individual at every locus)
  all.chrom.bd <- select(tidy.vcf, c(locus, CHROM, POS)) %>% distinct()
  names(all.chrom.bd) <- c("locus", "CHROM", "start")

  ## identify windows separately on each chromosome
  window.list <- lapply(unique(all.chrom.bd$CHROM), function(l){
    print(l) #printing each locus gives an indication of progress
    bd <- filter(all.chrom.bd, CHROM == l) %>% arrange(start) # make sure SNPs are in order by position

    snp.windows <- data.frame(do.call(bind_rows,lapply(1:nrow(bd), function(i){
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
  names(window.list) <- unique(all.chrom.bd$CHROM)
  
  return(window.list)
}
