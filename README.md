### Analysis can be found [HERE](https://github.com/Rajan-sust/covid19-Assembly/blob/development/script/Assembly_analysis.md).

### Scripts for assembly

```
script/paired_assembler.sh
```

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
conda install -c bioconda soapdenovo2
conda install -c bioconda minia
```

```
wget https://sourceforge.net/projects/quast/files/quast-5.0.2.tar.gz
tar -zxvf quast-5.0.2.tar.gz
```

### Required files for running the PE script
- [PE_561samples_final_561runs.txt](https://github.com/Rajan-sust/covid19-Assembly/blob/master/files/PE_561samples_final_561runs.txt)

- [adapter_pe.fa](https://github.com/Rajan-sust/covid19-Assembly/blob/master/files/adapter_pe.fa)

- [configSoapDenovo.txt](https://github.com/Rajan-sust/covid19-Assembly/blob/master/files/configSoapDenovo.txt)

- [MN908947.3.fasta](https://github.com/Rajan-sust/covid19-Assembly/blob/master/files/MN908947.3.fasta)

### Run the Script for PE

```
bash paired_assembler.sh PE_561samples_final_561runs.txt adapter_pe.fa configSoapDenovo.txt 
```

### SE


```
bash se_assembler.sh sra_sample.txt adapter_SE.fa
```
