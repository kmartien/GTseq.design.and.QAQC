#This is a step that you need to do for each population individually before running the Fst calculation in sliding windows
#It's done with angsd (which hopefully your cluster has, or you can install it with miniconda)
#!!!!!To run this first script you need to be located in the directory where your actual bam files are!!!!


#This is the first script we have to run 
#you can copy everything below here to submit through slum until the reach the next block of hashtags (where I wrote stop copying)
#!/bin/bash
#SBATCH -J saf_mden_pops      ## you can give whatever name you want here to identify your job
#SBATCH -o saf_mden_pops.log.out ## name a file to save log of your job
#SBATCH -e saf_mden_pops.error       ## save error message if any
#SBATCH --mail-user=morgan@bio.ku.dk  ## your email account to receive notification of job status
#SBATCH --mail-type=ALL 
#SBATCH -t 48:00:00    #Time to run it (hours, minutes, seconds)
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=8            # Number of CPU cores per task
#SBATCH --mem=20gb
############################################
#Load the module
module load angsd

#This is the reference genome you used for mapping (use the absolute path)
ref=/projects/mjolnir1/people/nrb613/mapped_bams/mden.reference.april.5/mMesDen1.pri.cur.20220330.fasta
#dir is where we will find the file containing a list of the bam files to be used. In my case I ran this where the bam files were but had the bam.list (containing the files to use) in a different directory (the calc_Fst directory). So the script will find the file that contains the names of bams to be used for a population in the directory indicated, but by default angsd will look for the actual files in the directory we are currently working in. Sorry that's confusing
dir=/projects/mjolnir1/people/nrb613/mapped_bams/raw.bams/high.quality.bam.files/angsd/Fst/mden/calc_Fst
#This is the output directory that is sort of the stem for all the analysis. It should be a directory where you run all the Fst analysis and contains bam.list files that have the bam files per location. I recommend keeping the calc_Fst ending because the script may refer back to it later on
out_dir=/projects/mjolnir1/people/nrb613/mapped_bams/raw.bams/high.quality.bam.files/angsd/Fst/mden/calc_Fst
#These are the chromosomes/scaffolds you are interested in. I have just isoloated the autosomes
#For me this file looks like this (just for the first 5) These should match the chromosome/scaffold headers in your reference genome
#SUPER_1
#SUPER_2
#SUPER_3
#SUPER_4
#SUPER_5
chrs=/projects/mjolnir1/people/nrb613/mapped_bams/mden.reference.april.5/mMesDen1.pri.cur.20220330.autosomes.txt

#Now we start the loop after we have defined some global variables
#The file unique.locations.txt looks like this
#Indo_Haw
#Indo_Afr
#Indo_Sou
#Atl_East
#Atl_Oth
#Atl_Bah

for pop in $(cat /projects/mjolnir1/people/nrb613/mapped_bams/raw.bams/high.quality.bam.files/angsd/Fst/mden/unique.locations.txt)
do
#Here num is grabbing the number of individuals in a population.
#For example the Indo_Haw.bam.list file will look like this
#Mde_L1-4_dedup_noRepeats.bam
#Mde_L2-9_dedup_noRepeats.bam
#Mde_L4-27_dedup_noRepeats.bam
#Mde_L4-29_dedup_noRepeats.bam
#Mde_L4-32_dedup_noRepeats.bam
#Mde_L5-34_dedup_noRepeats.bam
#So num is a variable that reads out the file (cat) then it counts the number of entries and cuts the output based on fields 6 and 7 usng a space delimeter " ". It's just counting the number individuals in a population so we can use half that number as a filter in the next step
num=`cat ${dir}/${pop}.bam.list | wc | cut -f6,7 -d " "`

#This halfnum is a variable that takes the number of individuals in a file (num) and divides it by two. This way angsd will only consider sites that had data for at least 50% of the individuals within a population
halfnum=$(expr ${num} / 2)

#Now angsd will look for each ${pop}.bam.list and run the -dosaf 1 command with it as an input. 
#The -doSaf 1 command creates SAF files that are files that contain sample allele frequencies
#We also include the following filters/options 
#-anc ${ref} (this places the reference genome that we mapped everything to in the command. It's previously been saved as a variable before the loop)
#-gl 1 (use the samtools algorithm for genotype likelihoods)
#-minMapQ 30 (use a minimum mapping score of 30)
#-minQ 30 (use a minimum base quality score of 30)
#-minInd ${halfnum} (only consider sites where we had data for at least half of the samples in the bam list)
#-rf ${chrs} (only consider the chromosomes indicated here. I have only used the autosomes and excluded the sex chromosomes)
#-nThreads 8 (this is the number of threads)
#-out ${out_dir}${pop}/${pop} (this will put the .saf output in a directory named after the population and the prefix of the saf file will also be the name of the population (e.g. Indo_Haw/Indo_Haw.saf.idx). Note you don't need to add the suffix ".saf.idx".  
angsd -bam ${dir}/${pop}.bam.list -anc ${ref} -dosaf 1 -gl 1 -minMapQ 30 -minQ 30 -minInd ${halfnum} -rf ${chrs} -nThreads 8 -out ${out_dir}/${pop}/${pop}
done

#Stop copying the first script just above this line
############################################
############################################
############################################
############################################
############################################
############################################
############################################
############################################

#This is the second script we have to run. At this point navigate to your calc_Fst directory 

#!/bin/bash
#SBATCH -J fst.mden      ## you can give whatever name you want here to identify your job
#SBATCH -o fst.mden.log.out ## name a file to save log of your job
#SBATCH -e fst.mden.error       ## save error message if any
#SBATCH --mail-user=morgan@bio.ku.dk  ## your email account to receive notification of job status
#SBATCH --mail-type=ALL 
#SBATCH -t 20:00:00    
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=24            # Number of CPU cores per task
#SBATCH --mem=20gb
############################################
#This analysis relies on angsd so make sure it's downloaded on your cluster (can also load it through miniconda)
module load angsd 

#Rather than having to run through all the population comparisons manually, I made a file that has also the population IDs that I want to comapre
#An example of how this looks is below (without the hashtags)
#Indo_Haw.Indo_Afr
#Indo_Haw.Indo_Sou
#Indo_Haw.Atl_East
#Indo_Haw.Atl_Oth
#Indo_Haw.Atl_Bah
#Indo_Afr.Indo_Sou
#Indo_Afr.Atl_East
#Indo_Afr.Atl_Oth
#Indo_Afr.Atl_Bah
#Indo_Sou.Atl_East
#Indo_Sou.Atl_Oth
#Indo_Sou.Atl_Bah
#Atl_East.Atl_Oth
#Atl_Oth.Atl_Bah

#So in this case I only wanted to run Fst in sliding windows between populations within the same ocean basin (not lose power including SNPs that seperate out ocean basins because these populations generally seperate no matter which sites I used - because the differentiation is so strong)
#For example I have the two populations Indo_Haw and Indo_Afr. In the pop.compare.txt file I seperate the two population names by the delimeter "." so in the file I have one line as Indo_Haw.Indo_Afr

#So now I'm starting a loop that will go over all the population comparisons I'm interested in (in the pop.compare.txt file)
for comparison in $(cat pop.compare.txt)
do
#Here I'm splling the two populations by their IDs using a ".". So pop1 will be a variable with the name of the first population and pop2 will be a variable with the name of the second population I'm interested in. In Indo_Haw.Indo_Afr example, pop1=Indo_Haw and pop2=Indo_Afr

#These variables will be needed because for each population comparison ANGSD (and the loop) will write in files to directories it creates based on the populations

pop1="$(echo $comparison | cut -f1 -d".")"
pop2="$(echo $comparison | cut -f2 -d".")"

#This estimates the 2d-SFS
realSFS -P 20 ${pop1}/${pop1}.saf.idx ${pop2}/${pop2}.saf.idx > paired_2dsfs/$comparison.2dsfs

#This step creates an index file which will be needed
realSFS fst index -whichFst 1 ${pop1}/${pop1}.saf.idx ${pop2}/${pop2}.saf.idx -P 20 -sfs paired_2dsfs/${comparison}.2dsfs -fstout binary/${comparison}

#This step calculates the fst between the two populations in general (it gives you two values based on different calculations)
realSFS fst stats binary/${comparison}.fst.idx -P 20 > final_fsts/${comparison}.fst

#This step calculates the fst in nonoverlapping windows of 100 bps
realSFS fst stats2 binary/${comparison}.fst.idx -win 100 -step 100 -P 20 > final_fsts/${comparison}.100bp.windows.fst

done