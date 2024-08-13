#!/bin/bash
#SBATCH -J genotype.likelihoods      ## you can give whatever name you want here to identify your job
#SBATCH -o genotype.likelihoods.log.out ## name a file to save log of your job
#SBATCH -e genotype.likelihoods.min.test.error       ## save error message if any
#SBATCH --mail-user=morgan@bio.ku.dk  ## your email account to receive notification of job status
#SBATCH --mail-type=ALL 
#SBATCH --time=12:00:00    
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=20            # Number of CPU cores per task
#SBATCH --mem=10gb
############################################
module load angsd

#This "bamlist" should contain a list of bams files and they should be located in the same directory that you're running this script from. OR they should contain the absolute path to the bam file.

#bamlist example
#EB101_dedup_noRepeats.bam
#EB102_dedup_noRepeats.bam
#EB103_dedup_noRepeats.bam
#EB105_dedup_noRepeats.bam
#EB106_dedup_noRepeats.bam
#EB107_dedup_noRepeats.bam
#EB108_dedup_noRepeats.bam
#EB109_dedup_noRepeats.bam

bamlist=/projects/mjolnir1/people/nrb613/X204SC22021885-Z01-F004/combined.nuclear/redo/E.b.bam.list.doc.1.txt

#This is just the directory where you want the beagle file to be written to. The part after the last "/" will be the name of the file (e.g. EB.genotype.likelihoods.doc.1.map30.autosomes.no.relatives) and the suffix will be appended to that file name.
out=/projects/mjolnir1/people/nrb613/X204SC22021885-Z01-F004/combined.nuclear/redo/EB.genotype.likelihoods.doc.1.autosomes/EB.genotype.likelihoods.doc.1.map30.autosomes.no.relatives

#This is a list of chromosomes/scaffolds that you want to use (I just isolated the autosomes)
chr=/projects/mjolnir1/people/nrb613/X204SC22021885-Z01-F004/reference.genome/Erignathus_barbatus_HiC_autosomes.txt
#e.g.
#HiC_scaffold_1
#HiC_scaffold_2
#HiC_scaffold_3
#HiC_scaffold_4
#HiC_scaffold_5

#Here inds is just the number of individuals you want to include in order for ANGSD to consider a site. So if you only want to call a SNP if there was data for at least 75% of the individuals than I would change inds to 72 if I had 96 individuals in total. Note that the number of individuals considered will be appended to the output file name so it's easier to keep track of if you multiple iterations of this. I usually try 95%, 75%, 50% to see if the results change. 

#The rest of the filters/options are
#-GL 1 use the samtools genotype likelihoods algorithm
#-rf is just the list of chromosomes to be used (in a text file)
#-doMaf 1 calculate major and minor allele frequencies
#-SNP_pval remove sites with a pvalue larger than 1e-6. This seems standard for the confidence that there is actually a SNP at a given position
#-minMaf 0.05 This removes SNPs that are super rare (minor allele frequency less than 5%)
#-doGlf 2 Print out the beagle likelihood file 
#-doMajorMinor 1 Infer the major and minor alleles from the likelihoods
#-minInd only consider a site if there is data present for ${inds} number of individuals (adjust this below). I think 75% of individuals is a good starting point. But remember to put the number of actual individuals in and not 0.75 (for 75%)
#-bam this is your bamlist (contains a list of all the mapped files you're interested in)
#-uniqueonly 1 Consider only reads that mapped to one location on the reference genome
#-only_proper_pairs 1 Ignore orphaned reads (that don't have their mate pair)
#-skipTriallelic 1 Only consider biallelic SNPs
#-minMapQ 30 Minimum read mapping score of 30
#-minQ Minimum base quality score of 30


#75%
inds=72
angsd -GL 1 -rf ${chr} -out ${out}.${inds} -doMaf 1 -SNP_pval 1e-6 -nThreads 20 -minMaf 0.05 -doGlf 2 -doMajorMinor 1 -minInd ${inds} -bam ${bamlist} -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -skipTriallelic 1 -minMapQ 30 -minQ 30