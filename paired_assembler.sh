#!/bin/bash
# $1 -> path of the paired_sra_list.txt
# $2 -> path of the adapter file for PE
# assemblers with multi k-mer assembly approach: spades, metaspades, megahit, trinity
# assemblers with fixed k-mer assembly approach: velvet, abyss 


mapfile -t sra_list < "${1}"

mkdir -p -v fastqc multiqc trimmomatic trimmed_fastqc trimmed_multiqc \
            megahit trinity velvet

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
  python quast-5.0.2/quast.py -o "quast_results/megahit/${sra}" -r MN908947.3.fasta -t 40 "megahit/${sra}/${sra}.contigs.fa"
  
  Trinity --seqType fq --max_memory 10G --left  "${pgz_1}" --right "${pgz_2}" --no_bowtie --CPU 40 --full_cleanup
  mv trinity_out_dir.Trinity.fasta "trinity/${sra}.fasta"
  python quast-5.0.2/quast.py -o "quast_results/trinity/${sra}" -r MN908947.3.fasta -t 40 "trinity/${sra}.fasta"

  mkdir -p "abyss/${sra}-63"
  cd "abyss/${sra}-63/"
  abyss-pe name="${sra}" k=63 j=40 in="../../${pgz_1} ../../${pgz_2}"
  cd ../..
  python quast-5.0.2/quast.py -o "quast_results/abyss/${sra}-63" -r MN908947.3.fasta -t 40 "abyss/${sra}-63/${sra}-scaffolds.fasta"


  mkdir -p "abyss/${sra}-127"
  cd "abyss/${sra}-127/"
  abyss-pe name="${sra}" k=127 j=40 in="../../${pgz_1} ../../${pgz_2}"
  cd ../..
  python quast-5.0.2/quast.py -o "quast_results/abyss/${sra}-127" -r MN908947.3.fasta -t 40 "abyss/${sra}-127/${sra}-scaffolds.fasta"

  mkdir -p "spades/${sra}"
  spades.py -1 "${pgz_1}" -2 "${pgz_2}" --rna -t 40 -o "spades/${sra}"
  python quast-5.0.2/quast.py -o "quast_results/spades/${sra}" -r MN908947.3.fasta -t 40 "spades/${sra}/transcripts.fasta"
  

  mkdir -p "metaspades/${sra}"
  metaspades.py -1 "${pgz_1}" -2 "${pgz_2}" -t 40 -o "metaspades/${sra}"
  python quast-5.0.2/quast.py -o "quast_results/metaspades/${sra}" -r MN908947.3.fasta -t 40 "metaspades/${sra}/scaffolds.fasta"
  

  velveth "velvet/${sra}" 31 -short -separate -fastq "${pgz_1}" "${pgz_2}"
  velvetg "velvet/${sra}" -read_trkg yes
  python quast-5.0.2/quast.py -o "quast_results/velvet/${sra}" -r MN908947.3.fasta -t 40 "velvet/${sra}/contigs.fa"

done

multiqc -f fastqc -o multiqc
multiqc -f trimmed_fastqc -o trimmed_multiqc