#!/bin/bash
# $1 -> path of the paired_sra_list.txt
# $2 -> path of the adapter file for PE
# velvet, abyss and soapdenovo2: kmer = 21,63, 99 
# quast will not create output if contifs are <500bp

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

# abyss21
  mkdir -p "abyss21/${sra}"
  cd "abyss21/${sra}/"
  abyss-pe name="${sra}" k=63 j=40 in="../../${pgz_1} ../../${pgz_2}"
  rm *dot*
  cd ../..
  quast_cp_zip_rm "abyss21" "${sra}" "abyss21/${sra}/${sra}-scaffolds.fa"

# abyss63
  mkdir -p "abyss63/${sra}"
  cd "abyss63/${sra}/"
  abyss-pe name="${sra}" k=63 j=40 in="../../${pgz_1} ../../${pgz_2}"
  rm *dot*
  cd ../..
  quast_cp_zip_rm "abyss63" "${sra}" "abyss63/${sra}/${sra}-scaffolds.fa"
  
# abyss99
  mkdir -p "abyss99/${sra}"
  cd "abyss99/${sra}/"
  abyss-pe name="${sra}" k=99 j=40 in="../../${pgz_1} ../../${pgz_2}"
  rm *dot*
  cd ../..
  quast_cp_zip_rm "abyss99" "${sra}" "abyss99/${sra}/${sra}-scaffolds.fa"  

# spades
  mkdir -p "spades/${sra}"
  spades.py -1 "${pgz_1}" -2 "${pgz_2}" --rna -t 40 -o "spades/${sra}"
  quast_cp_zip_rm "spades" "${sra}" "spades/${sra}/transcripts.fasta"
  
# metaspades
  mkdir -p "metaspades/${sra}"
  metaspades.py -1 "${pgz_1}" -2 "${pgz_2}" -t 40 -o "metaspades/${sra}"
  quast_cp_zip_rm "metaspades" "${sra}" "metaspades/${sra}/scaffolds.fasta"

#  mkdir - p velvet
#  velveth "velvet/${sra}" 21 -short -separate -fastq "${pgz_1}" "${pgz_2}"
#  velvetg "velvet/${sra}" -read_trkg yes
#  rm "velvet/${sra}/*Graph*"
#  quast_cp_zip_rm "velvet" "${sra}" "velvet/${sra}/contigs.fa"
  
#  velveth "metavelvet/${sra}" 21 -short -separate -fastq "${pgz_1}" "${pgz_2}"
#  velvetg "metavelvet/${sra}" -read_trkg yes
#  meta-velvetg "metavelvet/${sra}"
#  rm "metavelvet/${sra}/*Graph*"
#  quast_cp_zip_rm "metavelvet" "${sra}" "metavelvet/${sra}/contigs.fa"

# ray21
  mkdir -p "ray21"
  cd "ray21"
  gunzip ../"${pgz_1}" ../"${pgz_2}"
  mpiexec -n 10 Ray -k 21 -p  ../"trimmomatic/${sra}_1P.fastq" ../"trimmomatic/${sra}_2P.fastq" -o "${sra}"
  rm -r "${sra}"/BiologicalAbundances/
  cd ../
  gzip "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_2P.fastq"
  quast_cp_zip_rm "ray21" "${sra}" "ray21/${sra}/Scaffolds.fasta"

# ray63
  mkdir -p "ray63"
  cd "ray63"
  gunzip ../"${pgz_1}" ../"${pgz_2}"
  mpiexec -n 10 Ray -k 63 -p  ../"trimmomatic/${sra}_1P.fastq" ../"trimmomatic/${sra}_2P.fastq" -o "${sra}"
  rm -r "${sra}"/BiologicalAbundances/
  cd ../
  gzip "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_2P.fastq"
  quast_cp_zip_rm "ray63" "${sra}" "ray63/${sra}/Scaffolds.fasta"

# ray99
  mkdir -p "ray99"
  cd "ray99"
  gunzip ../"${pgz_1}" ../"${pgz_2}"
  mpiexec -n 10 Ray -k 99 -p  ../"trimmomatic/${sra}_1P.fastq" ../"trimmomatic/${sra}_2P.fastq" -o "${sra}"
  rm -r "${sra}"/BiologicalAbundances/
  cd ../
  gzip "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_2P.fastq"
  quast_cp_zip_rm "ray99" "${sra}" "ray99/${sra}/Scaffolds.fasta"


  rm -rf "data/"
  rm -rf "trimmomatic/"
  
done

metaquast.py "output/*_PE.fasta" -o "output/meta_quast_out" -r MN908947.3.fasta -g Sars_cov_2.ASM985889v3.100.gff3 -t 40 --silent
multiqc -f fastqc -o multiqc
multiqc -f trimmed_fastqc -o trimmed_multiqc