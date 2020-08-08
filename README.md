### Analysis can be found [HERE](https://github.com/Rajan-sust/covid19-Assembly/blob/development/script/Assembly_analysis.md).

### Single and Paired SRA list files

```
files/PE_561samples_final_561runs.txt
files/SE_300samples_final_300runs.txt
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