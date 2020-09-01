#!/bin/bash

mapfile -t sra_list < "${1}"

mkdir -p -v fastqc trimmed_fastqc trimmomatic

for sra in "${sra_list[@]}"; do
  echo "started running analysis for ${sra}"
  fastq-dump --split-files --outdir data --gzip "${sra}"
  fastq_1="data/${sra}_1.fastq.gz"
  
  # fastqc
  fastqc -o fastqc "${fastq_1}"

  # trimmomatic
  fsgz="trimmomatic/${sra}_1S.fastq.gz"
  trimmomatic SE "${fastq_1}" "${fsgz}" SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:"${2}":2:40:15
  
  # fastqc after trimmomatic
  fastqc -o trimmed_fastqc "${fsgz}"

  # megahit
  # output file: megahit/${sra}/${sra}.contigs.fa
  mkdir -p megahit
  megahit -r "${fsgz}" -o "megahit/${sra}/" --out-prefix "${sra}"
  rm -rf "megahit/${sra}/intermediate_contigs/"

  # spades
  # output file: spades/${sra}/transcripts.fasta
  mkdir -p "spades/${sra}/"
  spades.py -s "${fsgz}" -o "spades/${sra}/" --rna

  # trinity
  # output file: trinity/${sra}/${sra}.fasta
  mkdir -p "trinity/${sra}/"
  Trinity --seqType fq --max_memory 10G --single "${fsgz}" --no_bowtie --CPU 40 --full_cleanup
  mv trinity_out_dir.Trinity.fasta "trinity/${sra}/${sra}.fasta"
  
  # velvet
  # output file: velvet${kmer}/${sra}/contigs.fa
  for kmer in 21 63 99; do
    mkdir -p "velvet${kmer}"
    velveth "velvet${kmer}/${sra}/" "${kmer}" -short -fastq "${fsgz}"
    velvetg "velvet${kmer}/${sra}/" -read_trkg yes
  done

  # abyss
  # output file: abyss${kmer}/${sra}/${sra}-unitigs.fa
  for kmer in 21 63 99; do
    mkdir -p "abyss${kmer}/${sra}/"
    cd "abyss${kmer}/${sra}/"
    abyss-pe name="${sra}" k=$kmer j=40 se="../../${fsgz}"
    rm *dot*
    cd ../..
  done

  # minia
  # output filename: minia${kmer}/${sra}/${sra}_1S.fastq.contigs.fa
  for kmer in 21 63 99; do
    mkdir -p "minia${kmer}/${sra}/"
    cd "minia${kmer}/${sra}/"
    minia -in "../../${fsgz}" -kmer-size "${kmer}" -abundance-min 3
    cd ../..
  done

  # ray
  gunzip "${fsgz}"
  for kmer in 21 63 99; do
    mkdir -p "ray${kmer}"
    cd "ray${kmer}"
    mpiexec -n 10 Ray -k $kmer -s "../trimmomatic/${sra}_1S.fastq" -o "${sra}"
    rm -r "${sra}/BiologicalAbundances/"
    cd ../
  done
  gzip "trimmomatic/${sra}_1S.fastq"

done
