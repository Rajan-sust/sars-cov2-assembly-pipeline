#!/bin/bash
# $1 -> path of the paired_sra_list.txt
# $2 -> path of the adapter file for PE
# velvet, abyss, ray and soapdenovo2: kmer = 21,63,99 
# quast will not create output if contigs are <500bp

mapfile -t sra_list < "${1}"

mkdir -p -v fastqc multiqc trimmomatic trimmed_fastqc trimmed_multiqc output

# ${1} -> assembler name; ${2} -> sra; ${3} -> start time
function keep_time_and_space() {
  end=$(date +%s.%N)
  runtime=$( echo "${end} - ${3}" | bc -l )
  echo "${1},${2},${runtime}" >> time.csv
  echo -e "---end---" >> space.txt
}

for sra in "${sra_list[@]}"; do
  echo "started running analysis for "${sra}""
  fastq-dump --split-files --outdir data --gzip "${sra}"
  fastq_1="data/${sra}_1.fastq.gz"
  fastq_2="data/${sra}_2.fastq.gz"
  
  fastqc -o fastqc "${fastq_1}" "${fastq_2}"
  
  pgz_1="trimmomatic/${sra}_1P.fastq.gz"
  ugz_1="trimmomatic/${sra}_1U.fastq.gz"
  pgz_2="trimmomatic/${sra}_2P.fastq.gz"
  ugz_2="trimmomatic/${sra}_2U.fastq.gz"

  trimmomatic PE "${fastq_1}" "${fastq_2}" "${pgz_1}" "${ugz_1}" "${pgz_2}" "${ugz_2}" \
  SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:"${2}":2:40:15

  fastqc -o trimmed_fastqc "${pgz_1}" "${pgz_2}"

  # megahit  
  mkdir -p megahit
  echo -e "---start---\nmegahit,${sra}" >> space.txt
  start=$(date +%s.%N)
  megahit -1 "${pgz_1}" -2 "${pgz_2}" -o "megahit/${sra}" --out-prefix "${sra}"
  keep_time_and_space "megahit" "${sra}" "${start}"
  
  # trinity
  mkdir -p "trinity/${sra}"
  echo -e "---start---\ntrinity,${sra}" >> space.txt
  start=$(date +%s.%N)
  Trinity --seqType fq --max_memory 10G --left  "${pgz_1}" --right "${pgz_2}" --no_bowtie --CPU 40 --full_cleanup
  keep_time_and_space "trinity" "${sra}" "${start}"
  mv trinity_out_dir.Trinity.fasta "trinity/${sra}/${sra}.fasta"

  # spades
  mkdir -p "spades/${sra}"
  echo -e "---start---\nspades,${sra}" >> space.txt
  start=$(date +%s.%N)
  spades.py -1 "${pgz_1}" -2 "${pgz_2}" --rna -t 40 -o "spades/${sra}"
  keep_time_and_space "spades" "${sra}" "${start}"
  
  # metaspades
  mkdir -p "metaspades/${sra}"
  echo -e "---start---\nmetaspades,${sra}" >> space.txt
  start=$(date +%s.%N)
  metaspades.py -1 "${pgz_1}" -2 "${pgz_2}" -t 40 -o "metaspades/${sra}"
  keep_time_and_space "metaspades" "${sra}" "${start}"
  
  # abyss
  for kmer in 21 63 99; do
    mkdir -p "abyss${kmer}/${sra}"
    echo -e "---start---\nabyss${kmer},${sra}" >> space.txt
    start=$(date +%s.%N)
    cd "abyss${kmer}/${sra}/"
    abyss-pe name="${sra}" k=$kmer j=40 in="../../${pgz_1} ../../${pgz_2}"
    cd ../..
    keep_time_and_space "abyss${kmer}" "${sra}" "${start}"
  done 

  # velvet and metavelvet
  for kmer in 21 63 99; do
    mkdir - p "velvet${kmer}"
    echo -e "---start---\nvelvet${kmer},${sra}" >> space.txt
    start=$(date +%s.%N)
    velveth "velvet${kmer}/${sra}" $kmer -short -separate -fastq "${pgz_1}" "${pgz_2}"
    velvetg "velvet${kmer}/${sra}" -read_trkg yes
    keep_time_and_space "velvet${kmer}" "${sra}" "${start}"
  
    mkdir - p "metavelvet${kmer}"
    echo -e "---start---\nmetavelvet${kmer},${sra}" >> space.txt
    start=$(date +%s.%N)
    velveth "metavelvet${kmer}/${sra}" $kmer -short -separate -fastq "${pgz_1}" "${pgz_2}"
    velvetg "metavelvet${kmer}/${sra}" -read_trkg yes
    meta-velvetg "metavelvet${kmer}/${sra}"
    keep_time_and_space "metavelvet${kmer}" "${sra}" "${start}"
  done

  # ray
  for kmer in 21 63 99; do
    mkdir -p "ray${kmer}"
    echo -e "---start---\nray${kmer},${sra}" >> space.txt
    cd "ray${kmer}"
    gunzip ../"${pgz_1}" ../"${pgz_2}"
    start=$(date +%s.%N)
    mpiexec -n 10 Ray -k $kmer -p  ../"trimmomatic/${sra}_1P.fastq" ../"trimmomatic/${sra}_2P.fastq" -o "${sra}"
    cd ../
    keep_time_and_space "ray${kmer}" "${sra}" "${start}"
    gzip "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_2P.fastq"
  done

  # soapdenovo
  # $3 -> configSoapDenovo.txt
  # configSoapDenovo.txt must be end with newline.
  for kmer in 21 63 99; do
    mkdir -p "soapdenovo${kmer}/${sra}"
    echo -e "---start---\nsoapdenovo${kmer},${sra}" >> space.txt
    cd "soapdenovo${kmer}/${sra}"
    cat "../../${3}" > temp.txt
    echo "q1=../../trimmomatic/${sra}_1P.fastq" >> temp.txt
    echo "q2=../../trimmomatic/${sra}_2P.fastq" >> temp.txt
    start=$(date +%s.%N)
    SOAPdenovo-63mer all -s temp.txt -o "${sra}" -K "${kmer}"
    cd ../..
    keep_time_and_space "soapdenovo${kmer}" "${sra}" "${start}"

  done

  # minia
  for kmer in 21 63 99; do
    mkdir -p "minia${kmer}/${sra}/"
    echo -e "---start---\nminia${kmer},${sra}" >> space.txt
    cd "minia${kmer}/${sra}/"
    cat "../../trimmomatic/${sra}_1P.fastq" "../../trimmomatic/${sra}_2P.fastq" > "${sra}_M.fastq"
    start=$(date +%s.%N)
    minia -in "${sra}_M.fastq" -kmer-size "${kmer}" -abundance-min 3
    cd ../..
    keep_time_and_space "minia${kmer}" "${sra}" "${start}"
  done

  
done

