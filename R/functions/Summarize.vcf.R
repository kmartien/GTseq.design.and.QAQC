summarize.vcf <- function(vcf.dir, res.dir, fname, res.name){

  vcf <- paste0(vcf.dir,"/",fname)
  res <- paste0(res.dir,"/",res.name)
  print(res)
  
# depth indv/locus
  vcftools.idepth(vcf, res.fn = res)
  vcftools.siteDepth(vcf, res)
  vcftools.genoDepth(vcf, res)

# missing data indv/locus
  vcftools.imiss(vcf, res)
  vcftools.siteMiss(vcf, res)

# allele freq/indv freq buden
  vcftools.freqBurden(vcf, res)
  vcftools.freq2(vcf, res)
  vcftools.singletons(vcf, res)
  vcftools.012(vcf, res)

# heterozygosity per individual
  vcftools.het(vcf, res)

# SNP call quality
  vcftools.siteQual(vcf, res)

visualize.vcf.ind(res.dir, res.name)
visualize.vcf.loci(res.dir, res.name)

}