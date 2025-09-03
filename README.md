### Abstract
Coronavirus Disease 2019 (COVID-19), caused by severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2), has become a global pandemic following its initial emergence in China. SARS-CoV-2 has a positive-sense single-stranded RNA virus genome of around 30Kb. Using next-generation sequencing technologies, a large number of SARS-CoV-2 genomes are being sequenced at an unprecedented rate and being deposited in public repositories. For the de novo assembly of the SARS-CoV-2 genomes, a myriad of assemblers is being used, although their impact on the assembly quality has not been characterized for this virus. In this study, we aim to understand the variabilities on assembly qualities due to the choice of the assemblers.We performed 6648 de novo assemblies of 416 SARS-CoV-2 samples using eight different assemblers with different k-mer lengths. We used Illumina paired-end sequencing reads and compared the assembly quality of those assemblers. We showed that the choice of assembler plays a significant role in reconstructing the SARS-CoV-2 genome. Two metagenomic assemblers, e.g. MEGAHIT and metaSPAdes, performed better compared with others in most of the assembly quality metrics including, recovery of a larger fraction of the genome, constructing larger contigs and higher N50, NA50 values, etc. We showed that at least 09\% (259/2873) of the variants present in the assemblies between MEGAHIT and metaSPAdes are unique to one of the assembly methods.Our analyses indicate the critical role of assembly methods for assembling SARS-CoV-2 genome using short reads and their impact on variant characterization. This study could help guide future studies to determine the best-suited assembler for the de novo assembly of virus genomes.

### Analysis can be found [HERE](https://github.com/Rajan-sust/covid19-Assembly/blob/master/script/markdown/Assembly_analysis.md).

### SARS-CoV-2 Genome Assembly Pipeline for Single-End (SE) and Paired-End (PE) Data

```bash
# Create a new conda environment with Python 3.8
conda create -n sars-cov2 python=3.8

# Activate the environment
conda activate sars-cov2

# Add required channels for bioinformatics tools
conda config --env --add channels conda-forge
conda config --env --add channels bioconda

# Install essential packages
conda install -y sra-tools fastqc fastp megahit seqkit

# Download and install QUAST for assembly quality assessment
wget https://github.com/ablab/quast/releases/download/quast_5.3.0/quast-5.3.0.tar.gz
tar -xzf quast-5.3.0.tar.gz
cd quast-5.3.0
./setup.py install
cd ..
```

**Usage:**
```bash
bash main.sh -a <SRA_ACCESSION> [-o <OUTPUT_DIR>] [-t <THREADS>]
```

**Options:**
- `-a`    SRA accession number (**required**)
- `-o`    Output directory (default: current directory)
- `-t`    Number of threads (default: 4)
- `-h`    Show help message

**Example:**
```bash
bash main.sh -a SRR1234567 -o ./analysis -t 8
```

### Citation
```
@article{10.1093/bib/bbab102,
    author = {Islam, Rashedul and Raju, Rajan Saha and Tasnim, Nazia and Shihab, Istiak Hossain and Bhuiyan, Maruf Ahmed and Araf, Yusha and Islam, Tofazzal},
    title = {Choice of assemblers has a critical impact on de novo assembly of SARS-CoV-2 genome and characterizing variants},
    journal = {Briefings in Bioinformatics},
    volume = {22},
    number = {5},
    pages = {bbab102},
    year = {2021},
    month = {04},
    issn = {1477-4054},
    doi = {10.1093/bib/bbab102},
    url = {https://doi.org/10.1093/bib/bbab102},
    eprint = {https://academic.oup.com/bib/article-pdf/22/5/bbab102/40261256/bbab102.pdf},
}

```