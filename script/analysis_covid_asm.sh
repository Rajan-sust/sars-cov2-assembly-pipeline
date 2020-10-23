#path


# make window of 50bp
cd /Volumes/rony/drive/asm/covid19-Assembly/files/
bedtools makewindows -g MN908947.3.sizes.txt -w 50 >MN908947.3.50bpWindows.bed
bedtools makewindows -g MN908947.3.sizes.txt -w 100 >MN908947.3.100bpWindows.bed

# make bed file of alignment from paf file
cd /Volumes/rony/drive/asm/covid19-Assembly/files/VAR/
for f in *.srt.paf; do less $f | awk '{print $6 "\t" $8 "\t" $9}' | sort -k1,1 -k2,2n >$f.bed;done

#make matrix 
for f in *paf.bed; do echo $f; intersectBed -a /Volumes/rony/drive/asm/covid19-Assembly/files/MN908947.3.50bpWindows.bed -b $f -wao | awk '{print $1"_"$2"_"$3 "\t" $0}' | awk '!seen[$1]++' | awk '{print $8}' >$f.50bp_overlap;done

#print zero out file names
for f in *paf.bed.50bp_overlap; do wc -l $f;done | awk '$1 == 0 {print "rm" "\t" $2}' >zero_out_remove.sh

sh zero_out_remove.sh

#all assemblers
less /Volumes/rony/drive/asm/covid19-Assembly/files/MN908947.3.50bpWindows.bed | awk '{print $1"_"$2"_"$3}' >A.ID.all.50bp_overlap

((echo *megahit*.50bp_overlap |tr ' ' '\t') && (paste *megahit*.50bp_overlap)) >../megahit_alignment_matrix_50bp.tsv

#megahit
less /Volumes/rony/drive/asm/covid19-Assembly/files/MN908947.3.50bpWindows.bed | awk '{print $1"_"$2"_"$3}' >A.ID.megahit.50bp_overlap

((echo *megahit*.50bp_overlap |tr ' ' '\t') && (paste *megahit*.50bp_overlap)) >../megahit_alignment_matrix_50bp.tsv

#metaspades
less /Volumes/rony/drive/asm/covid19-Assembly/files/MN908947.3.50bpWindows.bed | awk '{print $1"_"$2"_"$3}' >A.ID.metaspade.50bp_overlap

((echo *metaspade*.50bp_overlap |tr ' ' '\t') && (paste *metaspade*.50bp_overlap)) >../metaspade_alignment_matrix_50bp.tsv

#trinity
less /Volumes/rony/drive/asm/covid19-Assembly/files/MN908947.3.50bpWindows.bed | awk '{print $1"_"$2"_"$3}' >A.ID.trinity.50bp_overlap

((echo *trinity*.50bp_overlap |tr ' ' '\t') && (paste *trinity*.50bp_overlap)) >../trinity_alignment_matrix_50bp.tsv

#abyss99
less /Volumes/rony/drive/asm/covid19-Assembly/files/MN908947.3.50bpWindows.bed | awk '{print $1"_"$2"_"$3}' >A.ID.abyss99.50bp_overlap

((echo *abyss99*.50bp_overlap |tr ' ' '\t') && (paste *abyss99*.50bp_overlap)) >../abyss99_alignment_matrix_50bp.tsv

#
less /Volumes/rony/drive/asm/covid19-Assembly/files/MN908947.3.50bpWindows.bed | awk '{print $1"_"$2"_"$3}' >A.ID.abyss63.50bp_overlap

((echo *abyss63*.50bp_overlap |tr ' ' '\t') && (paste *abyss63*.50bp_overlap)) >../abyss63_alignment_matrix_50bp.tsv

#remove tmp files
for f in *.50bp_overlap; do rm $f;done 

# move fastq data of libs without R2 file
cd /projects/epigenomics3/temp/rislam/assembly/rajan/asm_pe/fastqc/without_R2_SRR/
ls *zip | paste - | tr '_' '\t' |awk '{print $1}' | sort | uniq -c | awk '$1 == 2{print $2}' >list_SRR_wit_R1_R2.txt
#
while read line; do echo $line; 
mv $line* ../;
done <list_SRR_wit_R1_R2.txt

# make read qc matrix
cd /projects/epigenomics3/temp/rislam/assembly/rajan/asm_pe/fastqc
for f in *zip; do unzip $f;done

for f in *RR*/fastqc_data.txt; do echo $f; grep "Total Sequences\|Sequences flagged as poor quality\|Sequence length\|%GC" $f; done | paste - - - - -  | awk '{gsub("/fastqc_data.txt", ""); print }' >read_QC_matrix.txt

# run metaquast
cd /projects/epigenomics3/temp/rislam/assembly/rajan/asm_pe/

metaquast.py output/*_PE.fasta -o ./meta_quast_out_PE_all/ -r MN908947.3.fasta -g Sars_cov_2.ASM985889v3.100.gff3 -t 40 --fast --silent

#compress dir
tar -zcvf quast_results.tar.gz quast_results
#uncompress
tar xvzf myfolder.tar.gz

# bam files
for f in *bam; do samtools sort -o $f.sorted.bam $f;done
for f in *sorted.bam; do samtools index $f;done

## minimar2 data process
# sam to bam
cd /Users/rashedulislam/Documents/research/asm/bam/minimap2_alignment/sam_contigs/
for f in *sam; do echo $f; 
less $f | samtools view -bS | samtools sort -o $f.sorted.bam;
samtools index $f.sorted.bam;
done

mv *.sorted.bam* ../bam_contigs/

#picard dups marking
bam=/Users/rashedulislam/Documents/research/asm/bam/bwa_alignment/disagreement_MetaSpades_Megahit
cd /Users/rashedulislam/Documents/Tools/
for f in /Users/rashedulislam/Documents/research/asm/bam/bwa_alignment/best_worst_12bams/*.sorted.bam;
do ls $f;
#java -jar picard.jar MarkDuplicates I=$bam/$f O=$bam/$f.dupsMarked M=$bam/$f.marked_dup_metrics.txt
java -jar picard.jar MarkDuplicates I=$f O=$f.dupsMarked.bam M=$f.marked_dup_metrics.txt
done

#
for f in /projects/epigenomics3/temp/rislam/temp/bwa_alignment/*/*.sorted.bam;
do ll $f;
picard MarkDuplicates I=$f O=$f.sorted_dupsMarked.bam M=$f.marked_dup_metrics.txt;
done



java -jar picard.jar MarkDuplicates I=SRR12182136.bam.sorted.bam O=SRR12182136.bam.sorted.bam.sorted_dupsMarked.bam M=SRR12182136.bam.sorted.bam.marked_dup_metrics.txt



#mutation summary:
MEGA_META: https://github.com/istiakshihab/x-genomics-alignment-tester/blob/master/6k_MegaMeta_Assembler_FASTA_AND_PROCESSED/Final_Match_Missmatch.csv

META_TRINITY: https://github.com/istiakshihab/x-genomics-alignment-tester/blob/master/6k_TrinMegaMeta_Assembler_All_Miismatch/meta-trinity-match-mismatch-summary.csv

MEGA_TRINITY: https://github.com/istiakshihab/x-genomics-alignment-tester/blob/master/6k_TrinMegaMeta_Assembler_All_Miismatch/mega_trinity_filtered_call_summary.csv

## vcfs

MEGA_META_MERGED_VCF: https://github.com/istiakshihab/x-genomics-alignment-tester/tree/master/6k_MegaMeta_Assembler_FASTA_AND_PROCESSED/MERGED_VCF

META_TRINITY_MERGED_VCF: https://github.com/istiakshihab/x-genomics-alignment-tester/tree/master/6k_TrinMegaMeta_Assembler_All_Miismatch/TRINITY_BAM/META_TRINITY_MERGED_VCF

MEGA_TRINITY_MERGED_VCF: https://github.com/istiakshihab/x-genomics-alignment-tester/tree/master/6k_TrinMegaMeta_Assembler_All_Miismatch/TRINITY_BAM/MEGA_TRINITY_MERGED_VCF





