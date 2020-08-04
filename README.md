### Single and Paired SRA ID Selection

```
tail -n +2 COVID19_14.06.20_metadata_subsampled.csv | tr ',' '\t' | cut -f2,5 | grep "SINGLE" | cut -f1 > single_sra_list.txt 

tail -n +2 COVID19_14.06.20_metadata_subsampled.csv | tr ',' '\t' | cut -f2,5 | grep "PAIRED" | cut -f1 > paired_sra_list.txt  
```

### Requirements Installation

```
apt-get update
apt-get install sra-toolkit
apt-get install abyss
```

```
pip install multiqc
```

```
conda install -c bioconda fastqc
conda install -c bioconda trimmomatic
conda install -c bioconda megahit
conda install -c bioconda trinity
conda install -c bioconda spades
conda install -c bioconda velvet
conda install -c bioconda ray
```

```
wget https://sourceforge.net/projects/quast/files/quast-5.0.2.tar.gz
tar -zxvf quast-5.0.2.tar.gz
```

### Run the Script for PE

```
bash paired_assembler.sh sample_paired_sra_list.txt adapter_pe.fa
```