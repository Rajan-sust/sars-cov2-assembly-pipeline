#!/bin/bash
# $1 -> path of the paired_sra_list.txt
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
done