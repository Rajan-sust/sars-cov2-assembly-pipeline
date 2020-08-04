#!/bin/bash
# $1 -> path of the paired_sra_list.txt
# $2 -> path of the adapter file for PE

mapfile -t sra_list < "${1}"

mkdir -p -v fastqc multiqc trimmomatic trimmed_fastqc trimmed_multiqc \
            megahit trinity velvet output

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
  cp "megahit/${sra}/${sra}.contigs.fa" "output/${sra}_megahit_PE.fasta"
  cp "quast_results/megahit/${sra}/report.tsv" "output/${sra}_megahit_PE_quast.tsv"
   

  Trinity --seqType fq --max_memory 10G --left  "${pgz_1}" --right "${pgz_2}" --no_bowtie --CPU 40 --full_cleanup
  mv trinity_out_dir.Trinity.fasta "trinity/${sra}.fasta"
  python quast-5.0.2/quast.py -o "quast_results/trinity/${sra}" -r MN908947.3.fasta -t 40 "trinity/${sra}.fasta"
  cp "trinity/${sra}.fasta" "output/${sra}_trinity_PE.fasta"
  cp "quast_results/trinity/${sra}/report.tsv" "output/${sra}_trinity_PE_quast.tsv"

  mkdir -p "abyss/${sra}-63"
  cd "abyss/${sra}-63/"
  abyss-pe name="${sra}" k=63 j=40 in="../../${pgz_1} ../../${pgz_2}"
  cd ../..
  python quast-5.0.2/quast.py -o "quast_results/abyss/${sra}-63" -r MN908947.3.fasta -t 40 "abyss/${sra}-63/${sra}-scaffolds.fasta"
  cp "abyss/${sra}-63/${sra}-scaffolds.fasta" "output/${sra}_abyss63_PE.fasta"
  cp "quast_results/abyss/${sra}-63/report.tsv" "output/${sra}_abyss63_PE_quast.tsv"

  mkdir -p "abyss/${sra}-127"
  cd "abyss/${sra}-127/"
  abyss-pe name="${sra}" k=127 j=40 in="../../${pgz_1} ../../${pgz_2}"
  cd ../..
  python quast-5.0.2/quast.py -o "quast_results/abyss/${sra}-127" -r MN908947.3.fasta -t 40 "abyss/${sra}-127/${sra}-scaffolds.fasta"
  cp "abyss/${sra}-127/${sra}-scaffolds.fasta" "output/${sra}_abyss127_PE.fasta"
  cp "quast_results/abyss/${sra}-127/report.tsv" "output/${sra}_abyss127_PE_quast.tsv"

  mkdir -p "spades/${sra}"
  spades.py -1 "${pgz_1}" -2 "${pgz_2}" --rna -t 40 -o "spades/${sra}"
  python quast-5.0.2/quast.py -o "quast_results/spades/${sra}" -r MN908947.3.fasta -t 40 "spades/${sra}/transcripts.fasta"
  cp "spades/${sra}/transcripts.fasta" "output/${sra}_spades_PE.fasta"
  cp "quast_results/spades/${sra}/report.tsv" "output/${sra}_spades_PE_quast.tsv"

  mkdir -p "metaspades/${sra}"
  metaspades.py -1 "${pgz_1}" -2 "${pgz_2}" -t 40 -o "metaspades/${sra}"
  python quast-5.0.2/quast.py -o "quast_results/metaspades/${sra}" -r MN908947.3.fasta -t 40 "metaspades/${sra}/scaffolds.fasta"
  cp "metaspades/${sra}/scaffolds.fasta" "output/${sra}_metaspades_PE.fasta"
  cp "quast_results/metaspades/${sra}/report.tsv" "output/${sra}_metaspades_PE_quast.tsv"

  velveth "velvet/${sra}" 31 -short -separate -fastq "${pgz_1}" "${pgz_2}"
  velvetg "velvet/${sra}" -read_trkg yes
  python quast-5.0.2/quast.py -o "quast_results/velvet/${sra}" -r MN908947.3.fasta -t 40 "velvet/${sra}/contigs.fa"
  cp "velvet/${sra}/contigs.fa" "output/${sra}_velvet_PE.fasta"
  cp "quast_results/velvet/${sra}/report.tsv" "output/${sra}_metaspades_PE_quast.tsv"

  
  mkdir -p "ray/${sra}"
  gunzip "${pgz_1}" "${pgz_2}"
  mpiexec -n 10 Ray -k 21 -p  "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_1P.fastq" -o "ray/${sra}"
  gzip "trimmomatic/${sra}_1P.fastq" "trimmomatic/${sra}_1P.fastq"

done

multiqc -f fastqc -o multiqc
multiqc -f trimmed_fastqc -o trimmed_multiqc