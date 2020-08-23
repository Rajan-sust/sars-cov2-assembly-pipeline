#!/bin/bash
# $1 -> path of the paired_sra_list.txt
# $2 -> path of the adapter file for PE
# velvet, abyss, ray and soapdenovo2: kmer = 21,63,99 
# quast will not create output if contigs are <500bp

mapfile -t sra_list < "${1}"

mkdir -p -v fastqc multiqc trimmomatic trimmed_fastqc trimmed_multiqc output

# $1: assembly name, $2: sra, $3: assembly output full path
function quast_cp_zip_rm() {
  quast.py -o "quast_results/${1}/${2}" -r MN908947.3.fasta -g Sars_cov_2.ASM985889v3.100.gff3 -t 40 "${3}" --space-efficient --no-plots --silent
  cp "${3}" "output/${2}_${1}_PE.fasta"
  cp "quast_results/${1}/${2}/report.tsv" "output/${2}_${1}_PE_quast.tsv"

  if [[ -e "${1}_backup.zip" ]]; then
    zip -ur "${1}_backup.zip" "${1}/${2}"
  else
    zip -r "${1}_backup.zip" "${1}/${2}"
  fi
  
  rm -rf "${1}/${2}"

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
  megahit -1 "${pgz_1}" -2 "${pgz_2}" -o "megahit/${sra}" --out-prefix "${sra}"
  rm -rf "megahit/${sra}/intermediate_contigs"
  quast_cp_zip_rm "megahit" "${sra}" "megahit/${sra}/${sra}.contigs.fa"
  
# trinity
  mkdir -p "trinity/${sra}"
  Trinity --seqType fq --max_memory 10G --left  "${pgz_1}" --right "${pgz_2}" --no_bowtie --CPU 40 --full_cleanup
  mv trinity_out_dir.Trinity.fasta "trinity/${sra}/${sra}.fasta"
  quast_cp_zip_rm "trinity" "${sra}" "trinity/${sra}/${sra}.fasta"

# spades
  mkdir -p "spades/${sra}"
  spades.py -1 "${pgz_1}" -2 "${pgz_2}" --rna -t 40 -o "spades/${sra}"
  quast_cp_zip_rm "spades" "${sra}" "spades/${sra}/transcripts.fasta"
  
# metaspades
  mkdir -p "metaspades/${sra}"
  metaspades.py -1 "${pgz_1}" -2 "${pgz_2}" -t 40 -o "metaspades/${sra}"
  quast_cp_zip_rm "metaspades" "${sra}" "metaspades/${sra}/scaffolds.fasta"
  
# abyss
  for kmer in 21 63 99; do
    mkdir -p "abyss${kmer}/${sra}"
    cd "abyss${kmer}/${sra}/"
    abyss-pe name="${sra}" k=$kmer j=40 in="../../${pgz_1} ../../${pgz_2}"
    rm *dot*
    cd ../..
    quast_cp_zip_rm "abyss${kmer}" "${sra}" "abyss${kmer}/${sra}/${sra}-scaffolds.fa"
  done 

# velvet and metavelvet
  for kmer in 21 63 99; do
    mkdir - p "velvet${kmer}"
    velveth "velvet${kmer}/${sra}" $kmer -short -separate -fastq "${pgz_1}" "${pgz_2}"
    velvetg "velvet${kmer}/${sra}" -read_trkg yes
    rm -f velvet${kmer}/${sra}/*Graph*
    quast_cp_zip_rm "velvet${kmer}" "${sra}" "velvet${kmer}/${sra}/contigs.fa"
  
    mkdir - p "metavelvet${kmer}"
    velveth "metavelvet${kmer}/${sra}" $kmer -short -separate -fastq "${pgz_1}" "${pgz_2}"
    velvetg "metavelvet${kmer}/${sra}" -read_trkg yes
    meta-velvetg "metavelvet${kmer}/${sra}"
    rm -f metavelvet${kmer}/${sra}/*Graph*
    quast_cp_zip_rm "metavelvet${kmer}" "${sra}" "metavelvet${kmer}/${sra}/contigs.fa"
  done

# ray
  for kmer in 21 63 99; do
    mkdir -p "ray${kmer}"
    cd "ray${kmer}"
    gunzip ../"${pgz_1}" ../"${pgz_2}"
    mpiexec -n 10 Ray -k $kmer -p  ../"trimmomatic/${sra}_1P.fastq" ../"trimmomatic/${sra}_2P.fastq" -o "${sra}"
    rm -r "${sra}"/BiologicalAbundances/
    cd ../
    # gzip "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_2P.fastq"
    quast_cp_zip_rm "ray${kmer}" "${sra}" "ray${kmer}/${sra}/Scaffolds.fasta"
  done

  # soapdenovo
  # $3 -> configSoapDenovo.txt
  for kmer in 21 63; do
    mkdir -p "soapdenovo/${sra}-${kmer}"
    cd "soapdenovo/${sra}-${kmer}/"
    cat "../../${3}" > temp.txt
    echo "q1=../../trimmomatic/${sra}_1P.fastq" >> temp.txt
    echo "q2=../../trimmomatic/${sra}_2P.fastq" >> temp.txt
    SOAPdenovo-63mer all -s temp.txt -o "${sra}" -K "${kmer}"
    rm -f temp.txt
    cd ../..
    quast_cp_zip_rm "soapdenovo" "${sra}-${kmer}" "soapdenovo/${sra}-${kmer}/${sra}.scafSeq"
  done

  # minia
  mkdir -p "minia/${sra}/"
  cd "minia/${sra}/"
  cat "../../trimmomatic/${sra}_1P.fastq" "../../trimmomatic/${sra}_2P.fastq" > "${sra}_M.fastq"
  minia -in "${sra}_M.fastq" -kmer-size 31 -abundance-min 3
  rm -f "${sra}_M.fastq"
  cd ../..
  quast_cp_zip_rm "minia" "${sra}" "minia/${sra}/${sra}_M.contigs.fa"


  rm -f "data/*"
  rm -f "trimmomatic/*"
  
done

rm -rf megahit trinity spades metaspades abyss[0-9][0-9] ray[0-9][0-9] velvet[0-9][0-9] metavelvet[0-9][0-9]

metaquast.py output/*.fasta -o "output/meta_quast_out" -r MN908947.3.fasta -g Sars_cov_2.ASM985889v3.100.gff3 -t 40 --silent

# giving error: raise ImportError('A recent version of Python 3 is required.')
#multiqc -f fastqc -o multiqc
#multiqc -f trimmed_fastqc -o trimmed_multiqc