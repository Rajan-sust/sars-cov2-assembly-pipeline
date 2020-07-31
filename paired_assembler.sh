#!/bin/bash
# $1 -> path of the paired_sra_list.txt
# $2 -> path of the adapter file for PE
mapfile -t sra_list < "${1}" # file to array
for sra in "${sra_list[@]}"; do
  fastq-dump --split-files --outdir data --gzip "${sra}"
  fastq_1="data/${sra}_1.fastq.gz"
  fastq_2="data/${sra}_2.fastq.gz"
  
  if ! [ -d fastqc ]; then
    mkdir fastqc
  fi
  
  fastqc -o fastqc "${fastq_1}" "${fastq_2}"
  
  if ! [ -d multiqc ]; then
    mkdir multiqc
  fi
  
  multiqc -f fastqc -o multiqc
  # rm -rf fastqc
  if ! [ -d trimmomatic ]; then
    mkdir trimmomatic
  fi

  pgz_1="trimmomatic/${sra}_1P.fastq.gz"
  ugz_1="trimmomatic/${sra}_1U.fastq.gz"
  pgz_2="trimmomatic/${sra}_2P.fastq.gz"
  ugz_2="trimmomatic/${sra}_2U.fastq.gz"
  
  trimmomatic PE "${fastq_1}" "${fastq_2}" "${pgz_1}" "${ugz_1}" "${pgz_2}" "${ugz_2}" \
  SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:"${2}":2:40:15

done