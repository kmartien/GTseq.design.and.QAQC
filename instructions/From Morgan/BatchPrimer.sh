#The first step for me was to take the known SNP position (e.g SUPER_1 990875) and extract a region around the SNP. I extracted 129 bps before and after the SNP I was interested in.

#I had a list of SNP sites in a file called "mden.snps.txt". An example for a few are listed below:
SUPER_1 990875
SUPER_1 2451315
SUPER_1 3052764
SUPER_1 3082394
SUPER_1 3795124
SUPER_1 3795363

#I needed to prepare a bed file that would allow me to extract the region surrounding the SNP using bedtools
#This command takes the SNPs, and defines the scaffold name, region (129bps) before the SNP, and region (129bs) after the SNP. It also ensures the three columns are all tab delimited.
cat mden.snps.txt | awk '{print $1"\t"$2-129"\t"$2+129}' > mden.bed.txt

#As an example, this is what the mden.bed.txt file looked like
SUPER_1	990746	991004
SUPER_1	2451186	2451444
SUPER_1	3052635	3052893
SUPER_1	3082265	3082523
SUPER_1	3794995	3795253

#Now that I have a bed file I can extract the actual sequence corresponding to each given window. This is performed with bedtools. I need to load bedtools, have access to the (-fi) reference genome, (-bed) my bed file, (-fo) and the name of my output file.
module load bedtools 
bedtools getfasta -fi mMesDen1.pri.cur.20220330.fasta -bed mden.bed.txt -fo mden.fasta

#The output of this step will look like this
>SUPER_1:990746-991004
TCATTTTATTCTAAAAATTAAATACCCTAAAAAGAAGAAAAAAATTTTTAAATAAGCTTCAAAAGTATGGGATCAGGAGAACTAGAGGTTTTCGCCAGTTTCATTTTCAACTGTGTTTCTGCCTTCAAATTTTACAATAGTTTCATGATCTTCTTAATTGTTCACACGCATGTTTTCAATCACCACTTCACATTTCCTAGAATTCATTCCTTTGGCATGTGCTGTTTATGGCCTACTGTGATTTATTTATCTGATCCC
>SUPER_1:2451186-2451444
CAGAGGAGTGACTGCATCTGAACAGCCTTGTGAGATGGTCCCAAAGCTTAATTATCCACAGTCAGCTAACCAGCTACTGCAGCTTTTAGTTAAATTGCAGTGGCTATACTTCCTGCATGTTCTGGTTACGACCTACTTCAAGAAAAGGATAGTAAAAGAGAAATACAGCTGCTTCCATATACGGAAATGAGGTCCCCAGAGAGCATATCTGTATCAAAGATGATGATTTTCAATAACAGAATAGAAAGGACGACAATA

#For the next step, I need to get a list of the headers we need to design primers for. For me I used the IDs of the regions (i.e. the headers from the fasta file). Here extract all of the lines starting with ">" and use what follows the first ">".  
cat mden.fasta | grep ">" | cut -f2 -d ">" > snp.headers.txt

#The output (snp.headers.txt) will look like this
SUPER_1:990746-991004
SUPER_1:2451186-2451444
SUPER_1:3052635-3052893
SUPER_1:3082265-3082523

for header in 

#Primer3 requires a standard input file. For each SNP header I have a loop that creates the input file needed to run primer3. Note at the very bottome of this loop, I have snp.headers.txt that we previously generated
#####Loop begin
while read -r header
do 
#This creates an ID of the SNP header (e.g. SUPER_1:990746-991004)
seqid=`echo ${header}`
#This gets the actual sequence previously extracted with bedtools
seq=`grep -A 1 -F "$header" mden.fasta | tail -n 1`

#Primer3 setup
#Now I create the input file needed to run Primer3. Note that I created a directory "primers", where all of the setup files are being written to. The name of the file will be the region name followed by setup.txt (e.g. SUPER_10:10267246-10267504.setup.txt)
#Note the PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION (interval list, default empty) option means not to design primers that overlap this region. So starting at 120 and for 20 bps afterwards, no primers will overlap with this region.
#I also restricted the produce size to be between 90-130 bps with the option PRIMER_PRODUCT_SIZE_RANGE=90-130.

Middle oligos may not overlap any region specified by this tag.
The associated value must be a space-separated list of
echo "SEQUENCE_ID=${seqid}" > primers/${header}.setup.txt
echo "SEQUENCE_TEMPLATE=${seq}" >> primers/${header}.setup.txt
echo "PRIMER_PRODUCT_SIZE_RANGE=90-130" >> primers/${header}.setup.txt
echo "PRIMER_INTERNAL_OLIGO_EXCLUDED_REGION=120,20" >> primers/${header}.setup.txt
echo "PRIMER_EXPLAIN_FLAG=1" >> primers/${header}.setup.txt
echo "=" >> primers/${header}.setup.txt

#This is the part that actually runs primer3 with the input file we generated for each SNP/region. Note that it's accessing the setup file and writing the primer3 output to the directory "primers"
primer3_core primers/${header}.setup.txt > primers/${header}.primer3.output.txt
done < snp.headers.txt
#####Loop end

#At this point I've switched into the directory "primers" and run the following commands.
#There will be multiple primer options in the output, but I wanted to extract the "best" primer, which is the one with a 0 suffix. 
#Now getting the top primer for each set

#Create a summary file
touch mden.primer3_pair_0_summary.txt
#Add the headers from the output files to the summary file we just created
echo -e "SEQUENCE_ID\tLEFT\tRIGHT\tSEQUENCE_TEMPLATE\tPRIMER_PAIR_0_PRODUCT_SIZE" > mden.primer3_pair_0_summary.txt

#This loops through each output and gathers the information for the desired primer and puts it into the summary file (mden.primer3_pair_0_summary.txt)
for file in *primer3.output.txt
do
SEQUENCE_ID=`cat ${file} | grep "SEQUENCE_ID" | cut -f2 -d "="` 
SEQUENCE_TEMPLATE=`cat ${file} | grep "SEQUENCE_TEMPLATE" | cut -f2 -d "="` 
PRIMER_PAIR_0_PRODUCT_SIZE=`cat ${file} | grep "PRIMER_PAIR_0_PRODUCT_SIZE" | cut -f2 -d "="` 
LEFT=`cat ${file} | grep "PRIMER_LEFT_0_SEQUENCE" | cut -f2 -d "="` 
RIGHT=`cat ${file} | grep "PRIMER_RIGHT_0_SEQUENCE" | cut -f2 -d "="`
echo -e "${SEQUENCE_ID}\t${LEFT}\t${RIGHT}\t${SEQUENCE_TEMPLATE}\t${PRIMER_PAIR_0_PRODUCT_SIZE}" >> mden.primer3_pair_0_summary.txt
done