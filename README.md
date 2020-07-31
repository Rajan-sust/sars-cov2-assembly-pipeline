### Single and Paired SRA ID Selection

```
tail -n +2 COVID19_14.06.20_metadata_subsampled.csv | tr ',' '\t' | cut -f2,5 | grep "SINGLE" | cut -f1 > single_sra_list.txt 

tail -n +2 COVID19_14.06.20_metadata_subsampled.csv | tr ',' '\t' | cut -f2,5 | grep "PAIRED" | cut -f1 > paired_sra_list.txt  
```

### Requirements Installation

```
conda install -c bioconda fastqc
pip install multiqc
conda install -c bioconda trimmomatic
```

### Run the Script for PE

```
bash paired_assembler.sh sample_paired_sra_list.txt adapter_pe.fa
```