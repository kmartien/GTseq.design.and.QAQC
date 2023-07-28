basic.vcf.filtering <- function(vcf.dir = "vcf", fname, res, minDP=5, minQ=20, meanDP=15, mac=3, rm.monomorphic=TRUE){
  
  fname <- paste0(vcf.dir, "/", fname)
  res <- paste0(vcf.dir, "/", res)
  filter.res <- list()
  
  #2.1 Remove genotypes with minimum depth < minDP (default = 5)
  filter.res$minDP <- vcftools.minDP(vcf.fn = fname, res.fn = paste0(fname, ".minDP", minDP), minDP)
  fname <- paste0(fname, ".minDP", minDP)
  print("minDP")
  print(filter.res$minDP)
  
  #2.2 Remove loci with minimum quality < minQ (default = 20)
  #  Includes only loci with QUAL > this threshold (calculated across individuals).
  filter.res$minQ <- vcftools.minQ(vcf.fn = paste0(fname, ".recode"), res.fn = paste0(fname, ".minQ", minQ), minQ)
  fname <- paste0(fname, ".minQ", minQ)
  print("minQ")
  print(filter.res$minQ)
  
  #2.3 Remove loci with mean depth > meanDP (default = 15)
  # Includes only loci with mean read depth across all individuals greater than or equal to the "min-meanDP" value
  filter.res$meanDP <- vcftools.meanDP(vcf.fn = paste0(fname, ".recode"), res.fn = paste0(fname, ".meanDP", meanDP), meanDP)
  fname <- paste0(fname, ".meanDP", meanDP)
  print("meanDP")
  print(filter.res$meanDP)
  
  #2.4 Remove loci with minor allele count < 3*
  #  Include only sites with Minor Allele Count greater than or equal to the "--mac" value. Allele count is simply the number of times that allele appears over all individuals at that site.
  if (mac > 0){
    filter.res$mac <- vcftools.mac(vcf.fn = paste0(fname, ".recode"), res.fn = paste0(fname, ".mac", mac), mac)
    fname <- paste0(fname, ".mac", mac)
    print("mac")
    print(filter.res$mac)
  }
  
  #2.5 Remove monomorphic sites*
  if (rm.monomorphic == TRUE){
    filter.res$monomorphic <- vcftools.monomorphic(vcf.fn = paste0(fname, ".recode"), res.fn = paste0(fname, ".monomorphic"))
    fname <- paste0(fname, ".monomorphic")
    print("monomorphic")
    print(filter.res$monomorphic)
  }
  return(list(filter.res = filter.res, fname = fname))
}