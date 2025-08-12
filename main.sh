accession='SRR11828424'
mkdir -p data

# fastq-dump --split-3 --origfmt --outdir data --gzip "${accession}"
adapter="files/adapter_pe.fa"

fastq_1="data/${accession}_1.fastq.gz"
fastq_2="data/${accession}_2.fastq.gz"

mkdir -p trimmomatic

f_paired="trimmomatic/${accession}_1_trimmed.fastq.gz"
f_unpaired="trimmomatic/${accession}_1_unpaired.fastq.gz"
r_paired="trimmomatic/${accession}_2_trimmed.fastq.gz"
r_unpaired="trimmomatic/${accession}_2_unpaired.fastq.gz"


# Trimmomatic: A flexible trimmer for Illumina Sequence Data.

# trimmomatic PE -phred33 \
# "${fastq_1}" "${fastq_2}" \
# "${f_paired}" "${f_unpaired}" "${r_paired}" "${r_unpaired}" \
# ILLUMINACLIP:${adapter}:2:40:15 \
# LEADING:3 \
# TRAILING:3 \
# SLIDINGWINDOW:4:20 \
# MINLEN:25

# megahit  
mkdir -p megahit
megahit -1 "${f_paired}" -2 "${r_paired}" -o "megahit/${accession}" --out-prefix "${accession}" --k-list 21,29,39,59,79,99

# rm -rf "megahit/${accession}/intermediate_contigs"
# quast_cp_zip_rm "megahit" "${accession}" "megahit/${accession}/${accession}.contigs.fa"



# conda install bioconda::sra-tools
# conda install bioconda::trimmomatic
# conda install bioconda::megahit

