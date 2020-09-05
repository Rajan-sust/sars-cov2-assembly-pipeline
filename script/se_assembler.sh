#!/bin/bash

mapfile -t sra_list < "${1}"

mkdir -p -v fastqc trimmed_fastqc trimmomatic output


# $1: assembly name, $2: sra, $3: assembly output full path
function quast_cp_zip_rm() {
  python quast-5.0.2/quast.py -o "quast_results/${1}/${2}/" -r MN908947.3.fasta -g Sars_cov_2.ASM985889v3.100.gff3 \
  -t 40 "${3}" --space-efficient --no-plots --silent
  
  cp "${3}" "output/${2}_${1}_SE.fasta"
  cp "quast_results/${1}/${2}/report.tsv" "output/${2}_${1}_SE_quast.tsv"

  if [[ -e "${1}_backup.zip" ]]; then
    zip -ur "${1}_backup.zip" "${1}/${2}/"
  else
    zip -r "${1}_backup.zip" "${1}/${2}/"
  fi
  
  # rm -rf "${1}/${2}/"

}

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
  quast_cp_zip_rm "megahit" "${sra}" "megahit/${sra}/${sra}.contigs.fa"

  # spades
  # output file: spades/${sra}/transcripts.fasta
  mkdir -p "spades/${sra}/"
  spades.py -s "${fsgz}" -o "spades/${sra}/" --rna
  quast_cp_zip_rm "spades" "${sra}" "spades/${sra}/transcripts.fasta"
  

  # metaspades
  # mkdir -p "metaspades/${sra}/"
  # metaspades.py -s "${fsgz}" -o "metaspades/${sra}/" -t 40

  # trinity
  # output file: trinity/${sra}/${sra}.fasta
  mkdir -p "trinity/${sra}/"
  Trinity --seqType fq --max_memory 10G --single "${fsgz}" --no_bowtie --CPU 40 --full_cleanup
  mv trinity_out_dir.Trinity.fasta "trinity/${sra}/${sra}.fasta"
  quast_cp_zip_rm "trinity" "${sra}" "trinity/${sra}/${sra}.fasta"

  # velvet
  # output file: velvet${kmer}/${sra}/contigs.fa
  for kmer in 21 63 99; do
    mkdir -p "velvet${kmer}"
    velveth "velvet${kmer}/${sra}/" "${kmer}" -short -fastq "${fsgz}"
    velvetg "velvet${kmer}/${sra}/" -read_trkg yes
    quast_cp_zip_rm "velvet${kmer}" "${sra}" "velvet${kmer}/${sra}/contigs.fa"
    
    # mkdir -p "metavelvet${kmer}"
    # velveth "metavelvet${kmer}/${sra}/" "${kmer}" -short -fastq "${fsgz}"
    # velvetg "metavelvet${kmer}/${sra}/" -read_trkg yes
    # meta-velvetg "metavelvet${kmer}/${sra}/"
    # rm -f "metavelvet${kmer}/${sra}/*Graph*"
  done

  # abyss
  # output file: abyss${kmer}/${sra}/${sra}-unitigs.fa
  for kmer in 21 63 99; do
    mkdir -p "abyss${kmer}/${sra}/"
    cd "abyss${kmer}/${sra}/"
    abyss-pe name="${sra}" k=$kmer j=40 se="../../${fsgz}"
    rm *dot*
    cd ../..
    quast_cp_zip_rm "abyss${kmer}" "${sra}" "abyss${kmer}/${sra}/${sra}-unitigs.fa"
  done

  # minia
  # output filename: minia${kmer}/${sra}/${sra}_1S.fastq.contigs.fa
  for kmer in 21 63 99; do
    mkdir -p "minia${kmer}/${sra}/"
    cd "minia${kmer}/${sra}/"
    minia -in "../../${fsgz}" -kmer-size "${kmer}" -abundance-min 3
    cd ../..
    quast_cp_zip_rm "minia${kmer}" "${sra}" "minia${kmer}/${sra}_1S.fastq.contigs.fa"
  done

  # ray
  gunzip "${fsgz}"
  for kmer in 21 63 99; do
    mkdir -p "ray${kmer}"
    cd "ray${kmer}"
    mpiexec -n 10 Ray -k $kmer -s "../trimmomatic/${sra}_1S.fastq" -o "${sra}"
    rm -rf "${sra}/BiologicalAbundances/"
    cd ../
    quast_cp_zip_rm "ray${kmer}" "${sra}" "ray${kmer}/${sra}/Scaffolds.fasta"
  done
  gzip "trimmomatic/${sra}_1S.fastq"
  
done