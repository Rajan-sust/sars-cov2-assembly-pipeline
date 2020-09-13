# make window of 50bp
cd ~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/
bedtools makewindows -g MN908947.3.sizes.txt -w 50 >MN908947.3.50bpWindows.bed

# make bed file of alignment from paf file
cd ~/Downloads/x-genomics-alignment-tester-master/tsv_excluded/metaspades/VAR/
for f in *.srt.paf; do less $f | awk '{print $6 "\t" $8 "\t" $9}' >$f.bed;done

#make matrix 
for f in *paf.bed; do echo $f; bedtools intersect -a ~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/MN908947.3.50bpWindows.bed -b $f -wao | awk '{print $1"_"$2"_"$3 "\t" $0}' | awk '!seen[$1]++' | awk '{print $8}' >$f.50bp_overlap;done

less ~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/MN908947.3.50bpWindows.bed | awk '{print $1"_"$2"_"$3}' >A.ID.50bp_overlap

((echo *.50bp_overlap |tr ' ' '\t') && (paste *.50bp_overlap)) >alignment_matrix_50bp.tsv

# run metaquast
cd /projects/epigenomics3/temp/rislam/assembly/rajan/asm_pe/

metaquast.py output/*_PE.fasta -o output/meta_quast_out/ -r MN908947.3.fasta -g Sars_cov_2.ASM985889v3.100.gff3 -t 40

#compress dir
tar -zcvf quast_results.tar.gz quast_results
#uncompress
tar xvzf myfolder.tar.gz

# 