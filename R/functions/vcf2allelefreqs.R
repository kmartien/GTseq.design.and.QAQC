vcf2allelefreqs <- function(tidy.vcf){
  tidy.alleles.long <- 
    select(tidy.vcf, c(CHROM, POS, locus, gt_GT)) %>% #select relevant columns
    mutate(allele = gsub("|", "/", gt_GT, fixed = TRUE)) %>% #standardize allele separator in genotypes
    separate_wider_delim(cols = "allele", delim = "/", names_sep = ".") %>% #separate genotypes ('C/T') into separate alleles ('C' and 'T')
    select(c(allele.1, allele.2, locus)) %>% #select relevant columns
    pivot_longer(cols = c(allele.1, allele.2), names_to = "allele") #reformat
  
  ## use table() to calculate frequencies at each locus 
  af2 <- do.call(bind_rows, lapply(unique(tidy.alleles.long$locus), function(l){
    filter(tidy.alleles.long, locus == l) %>% select(value) %>% table()
  }))
  af <- bind_cols(unique(tidy.alleles.long$locus), af2)
  names(af) <- c("locus", "allele.1", "allele.2")
  return(af)
}