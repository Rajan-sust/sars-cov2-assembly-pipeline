#!/bin/bash
# $1 -> path of the paired_sra_list.txt
# $2 -> path of the adapter file for PE
# assemblers with multi k-mer assembly approach: spades, metaspades, megahit, trinity
# assemblers with fixed k-mer assembly approach: velvet, abyss 

mapfile -t sra_list < "${1}"

mkdir -p -v fastqc multiqc trimmomatic trimmed_fastqc trimmed_multiqc \
            megahit trinity velvet output

# $1: assembly name, $2: sra, $3: assembly output full path
function quast_cp_zip_rm() {
  python quast-5.0.2/quast.py -o "quast_results/${1}/${2}" -r MN908947.3.fasta -t 40 "${3}"
  cp "${3}" "output/${2}_${1}_PE.fasta"
  cp "quast_results/${1}/${2}/report.tsv" "output/${2}_${1}_PE_quast.tsv"
  
  if [[ "${1}" == "trinity" ]]; then
    continue
  fi

  if [[ -e "${1}.zip" ]]; then
    zip -ur "${1}.zip" "${1}/${2}"
  else
    zip -r "${1}.zip" "${1}/${2}"
  fi
  
  rm -rf "${1}/${2}"

}

for sra in "${sra_list[@]}"; do
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

  
  megahit -1 "${pgz_1}" -2 "${pgz_2}" -o "megahit/${sra}" --out-prefix "${sra}"
  rm -rf "megahit/${sra}/intermediate_contigs"
  quast_cp_zip_rm "megahit" "${sra}" "megahit/${sra}/${sra}.contigs.fa"
  

  Trinity --seqType fq --max_memory 10G --left  "${pgz_1}" --right "${pgz_2}" --no_bowtie --CPU 40 --full_cleanup
  mv trinity_out_dir.Trinity.fasta "trinity/${sra}.fasta"
  quast_cp_zip_rm "trinity" "${sra}" "trinity/${sra}.fasta"
  
  mkdir -p "abyss63/${sra}"
  cd "abyss63/${sra}/"
  abyss-pe name="${sra}" k=63 j=40 in="../../${pgz_1} ../../${pgz_2}"
  cd ../..
  quast_cp_zip_rm "abyss63" "${sra}" "abyss63/${sra}/${sra}-scaffolds.fa"
  
  mkdir -p "abyss127/${sra}"
  cd "abyss127/${sra}/"
  abyss-pe name="${sra}" k=127 j=40 in="../../${pgz_1} ../../${pgz_2}"
  cd ../..
  quast_cp_zip_rm "abyss127" "${sra}" "abyss127/${sra}/${sra}-scaffolds.fa"
  
  mkdir -p "spades/${sra}"
  spades.py -1 "${pgz_1}" -2 "${pgz_2}" --rna -t 40 -o "spades/${sra}"
  quast_cp_zip_rm "spades" "${sra}" "spades/${sra}/transcripts.fasta"
  
  mkdir -p "metaspades/${sra}"
  metaspades.py -1 "${pgz_1}" -2 "${pgz_2}" -t 40 -o "metaspades/${sra}"
  quast_cp_zip_rm "metaspades" "${sra}" "metaspades/${sra}/scaffolds.fasta"
  
  velveth "velvet/${sra}" 31 -short -separate -fastq "${pgz_1}" "${pgz_2}"
  velvetg "velvet/${sra}" -read_trkg yes
  quast_cp_zip_rm "velvet" "${sra}" "velvet/${sra}/contigs.fa"
  
  mkdir -p "ray/${sra}"
  gunzip "${pgz_1}" "${pgz_2}"
  mpiexec -n 10 Ray -k 21 -p  "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_1P.fastq" -o "ray/${sra}"
  gzip "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_1P.fastq"
  quast_cp_zip_rm "ray" "${sra}" "ray/${sra}/Scafolds.fasta"

  # rm -f "data/${sra}*"
  
done

multiqc -f fastqc -o multiqc
multiqc -f trimmed_fastqc -o trimmed_multiqc