#!/bin/bash
clear
num_cores=4


rm -rf ./alignment
mkdir -p ./alignment


ids=$((for f in $(find ./fastq/by_id -name '*.fastq.gz'); do
  basename ${f} | grep -o '^z[[:digit:]]*'
done) | sort | uniq)


for i in ${ids}; do
  date
  echo --- Align ${i} ---
  
  bwa mem \
    -M \
    -t ${num_cores} \
    -R "@RG\tID:${i}\tLB:amplicon\tPL:ILLUMINA\tSM:${i}" \
    -v 3 \
    ./reference/fin_ref \
    $(find ./fastq/by_id -name "${i}*") | \
  samtools sort -@ ${num_cores} -O bam -o ./alignment/${i}.bam -
  echo " " 
done


echo --- Merge and Index BAM ---
samtools merge \
  -f \
  -o ./alignment/merged.bam \
  -@ ${num_cores} \
  $(find ./alignment -name "z*.bam")
samtools index \
  -@ ${num_cores} \
  ./alignment/merged.bam


echo --- Summarize Alignments ---
(
  (echo Number of alignments:; samtools view -c -@ ${num_cores} ./alignment/merged.bam) | paste - -;
  (echo Number of mapped reads:; samtools view -c -F 4 -@ ${num_cores} ./alignment/merged.bam) | paste - -;
  (echo Number of unmapped reads:; samtools view -c -f 4 -@ ${num_cores} ./alignment/merged.bam) | paste - -;
  (echo Number of mapped paired-end reads:; samtools view -c -f 1 -F 12 -@ ${num_cores} ./alignment/merged.bam) | paste - -
) > ./alignment/smry_alignment.txt


echo " "
date