---
title: "COV_assembly"
author: "Rashedul"
date: "7/14/2020"
output: 
  html_document: 
    keep_md: yes
---



### variant heatmap


```r
library(pheatmap)
library(tidyverse)

meta = read_csv("/Volumes/rony/drive/asm/covid19-Assembly/files/PE_561samples_final.csv")
mh = read.table("/Volumes/rony/drive/asm/covid19-Assembly/files/megahit_alignment_matrix_50bp.tsv", head = F)
ms = read.table("/Volumes/rony/drive/asm/covid19-Assembly/files/metaspade_alignment_matrix_50bp.tsv", head = F)
hs = cbind(mh, ms)
hst = data.frame(t(hs))
hst2 = hst %>% filter(X1 != "A.ID.megahit.50bp_overlap") %>% filter(X1 != "A.ID.metaspade.50bp_overlap")

hst3 = data.frame(str_split_fixed(hst2$X1, "_", 3), hst2)

megahit = hst3 %>% filter(hst3$X2.1 == "megahit") %>% select(-X2.1) %>% select(-X3.1) %>% select(-X1)
metaspades = hst3 %>% filter(hst3$X2.1 == "metaspades") %>% select(-X2.1) %>% select(-X3.1) %>% select(-X1)

both = merge(megahit, metaspades, by.x = "X1.1", by.y = "X1.1") 
both2 = inner_join(meta, both, by = c( "Run" = "X1.1"))
both3 = both2[, c(2,6, 15:ncol(both2))] 

#write.table(both3, "/Volumes/rony/drive/asm/covid19-Assembly/files/megahit_metaspades_heatmap_matrix.tsv", sep = "\t", quote = F, row.names = F, col.names = F)

both4 = data.frame(both3[,3:ncol(both3)]) 

#conver factor to numeric
both5 = both4 %>% 
  mutate_all(~as.numeric(as.character(.)))
rownames(both5) = both3$Run

#pheatmap(both5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, color = colorRampPalette(c("red", "white", "gray"))(1000), border_color = NA, gaps_col = 599)

#add annotation
anno = data.frame(Assay = both3$Assay_Type)
rownames(anno) = rownames(both5)

pheatmap(both5, annotation_row = anno, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, color = colorRampPalette(c("red", "white", "gray"))(1000), border_color = NA, gaps_col = 599, gaps_row = c(82, 82+65, 82+65+58, 82+65+58+88))
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

```r
plot.new()
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-1-2.png)<!-- -->

```r
pdf("/Volumes/rony/drive/asm/covid19-Assembly/plots/heatmap.pdf")
pheatmap(both5, annotation_row = anno, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, color = colorRampPalette(c("red", "white", "gray"))(1000), border_color = NA, gaps_col = 599, gaps_row = c(82, 82+65, 82+65+58, 82+65+58+88))
dev.off()
```

```
## pdf 
##   3
```

### assembly quality


```r
library(tidyverse)
library(reshape2)
library(viridis)
library(ggsci)

x = read_tsv("/Volumes/rony/drive/asm/covid19-Assembly/files/backup/bioRxiv_1074_assembly_report_PE_amplicon_viralRNA.tsv")
x = read_tsv("/Volumes/rony/drive/asm/covid19-Assembly/files/report_metaquast_PE_all.tsv")

x2 = melt(x, id.vars = "Assembly")
x2$Assembly = gsub("# ", "", x2$Assembly)
x3 = data.frame(str_split_fixed(x2$variable, "_", 3), value = as.numeric(x2$value), x2)

meta = read_csv("/Volumes/rony/drive/asm/covid19-Assembly/files/PE_561samples_final.csv")
meta2 = meta %>% select(Run, Assay_Type)
    
xm = left_join(x3, meta2, by = c("X1" = "Run"))

#
setwd("/Volumes/rony/drive/asm/covid19-Assembly/plots/")
make_boxplot = function(variableToPlot)
{
  x4 = xm %>% filter(Assembly == variableToPlot) 
  stat = x4 %>%
    group_by(X2) %>%
    summarize(mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE))
  
  print(stat)
  plot = ggplot(x4, aes(X2, value)) +
    geom_boxplot(aes(color = X2)) +
    geom_jitter(alpha = .5, width = .1, aes(color = X2)) +
    #scale_color_viridis(discrete=TRUE) +
    #scale_fill_material("red") +
    xlab("") +
    ylab(variableToPlot) +
    facet_grid(Assay_Type~.) +
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = .5),
        axis.text = element_text(color = "black", angle = 90, hjust = 1),
        strip.background =element_rect(fill="white"),
        legend.position = "none") 
  print(plot)
  
  filename = str_replace_all(variableToPlot, "[[:punct:]]", "")
  ggsave(filename=paste(filename,".pdf", sep="_"), width = 20, height = 30, units = "cm", device = 'pdf')
}

make_boxplot("Genome fraction (%)")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21       78.9  95.2 
##  2 abyss63       80.7  97.7 
##  3 abyss99       73.6  91.1 
##  4 megahit       92.2  99.7 
##  5 metaspades    91.8  99.7 
##  6 metavelvet21  17.9   8.49
##  7 metavelvet63  28.0  17.3 
##  8 metavelvet99  38.4  27.0 
##  9 ray21         87.0  98.1 
## 10 ray63         88.3  98.0 
## 11 ray99         88.3  98.0 
## 12 spades        92.4  99.7 
## 13 trinity       85.5  99.5 
## 14 velvet21      17.9   8.49
## 15 velvet63      28.2  17.3 
## 16 velvet99      38.8  26.3
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
make_boxplot("Largest contig")
```

```
## # A tibble: 16 x 3
##    X2             mean median
##    <fct>         <dbl>  <dbl>
##  1 abyss21      12037.  9464.
##  2 abyss63      11186.  6987 
##  3 abyss99       9199.  3306 
##  4 megahit      22662. 29840.
##  5 metaspades   21777. 27419 
##  6 metavelvet21  1015.   842 
##  7 metavelvet63  1520.  1091 
##  8 metavelvet99  2073.  1277 
##  9 ray21        18993. 20489 
## 10 ray63        18347. 19246 
## 11 ray99        18339. 19246 
## 12 spades       20494. 20681 
## 13 trinity      17099. 16048 
## 14 velvet21      1012.   838.
## 15 velvet63      1524.  1096 
## 16 velvet99      2086.  1261
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-2.png)<!-- -->

```r
make_boxplot("Total length")
```

```
## # A tibble: 16 x 3
##    X2              mean median
##    <fct>          <dbl>  <dbl>
##  1 abyss21       94283. 29022.
##  2 abyss63       53096. 29601 
##  3 abyss99       26531. 27970 
##  4 megahit      186413. 32970.
##  5 metaspades   199392. 32708.
##  6 metavelvet21  48649.  5008 
##  7 metavelvet63  46495.  7847 
##  8 metavelvet99  20185.  9235 
##  9 ray21         96334. 32158 
## 10 ray63         81436. 31157 
## 11 ray99         81434. 31172 
## 12 spades       337027. 37088.
## 13 trinity      207714. 40652.
## 14 velvet21      48826.  4826.
## 15 velvet63      46937.  8266 
## 16 velvet99      20247.  9397
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-3.png)<!-- -->

```r
make_boxplot("contigs")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21       92.5      6
##  2 abyss63       47.4      9
##  3 abyss99       18.3     10
##  4 megahit      160.      10
##  5 metaspades   178.       9
##  6 metavelvet21  76.1      8
##  7 metavelvet63  61.3     11
##  8 metavelvet99  24.1     12
##  9 ray21         75.0      9
## 10 ray63         63.6      9
## 11 ray99         63.7      9
## 12 spades       200.       9
## 13 trinity      199.      15
## 14 velvet21      76.4      8
## 15 velvet63      61.9     12
## 16 velvet99      24.1     12
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-4.png)<!-- -->

```r
make_boxplot("contigs (>= 1000 bp)")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21      20.5       3
##  2 abyss63      10.4       4
##  3 abyss99       5.38      3
##  4 megahit      40.3       4
##  5 metaspades   44.1       4
##  6 metavelvet21  2.99      0
##  7 metavelvet63  8.42      1
##  8 metavelvet99  4.10      1
##  9 ray21        22.9       5
## 10 ray63        18.5       5
## 11 ray99        18.5       5
## 12 spades       62.5       5
## 13 trinity      47.2       6
## 14 velvet21      2.97      0
## 15 velvet63      8.53      1
## 16 velvet99      4.16      1
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-5.png)<!-- -->

```r
make_boxplot("Largest alignment")
```

```
## # A tibble: 16 x 3
##    X2             mean median
##    <fct>         <dbl>  <dbl>
##  1 abyss21      14435. 15035 
##  2 abyss63      13145.  9994 
##  3 abyss99      10500.  4417 
##  4 megahit      23452. 29819 
##  5 metaspades   22533. 27419 
##  6 metavelvet21  1014.   776.
##  7 metavelvet63  1517.   970.
##  8 metavelvet99  2116.  1061 
##  9 ray21        18758. 21750 
## 10 ray63        18984. 21417 
## 11 ray99        18975. 21417 
## 12 spades       19660. 21999 
## 13 trinity      17494. 18538 
## 14 velvet21      1020.   792.
## 15 velvet63      1511.   982 
## 16 velvet99      2131.  1091
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-6.png)<!-- -->

```r
make_boxplot("mismatches per 100 kbp")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21       27.2   24.9
##  2 abyss63       26.8   24.6
##  3 abyss99       26.8   25.4
##  4 megahit       33.8   26.8
##  5 metaspades    33.1   26.8
##  6 metavelvet21  45.4   17.6
##  7 metavelvet63  33.2   22.2
##  8 metavelvet99  43.3   28.9
##  9 ray21         38.3   33.5
## 10 ray63         36.1   33.0
## 11 ray99         36.1   33.0
## 12 spades       198.    33.5
## 13 trinity       48.6   33.6
## 14 velvet21      44.5   19.0
## 15 velvet63      33.6   22.3
## 16 velvet99      42.6   26.9
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-7.png)<!-- -->

```r
make_boxplot("indels per 100 kbp")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21       6.30   3.38
##  2 abyss63       5.4    0   
##  3 abyss99       3.80   0   
##  4 megahit       4.00   0   
##  5 metaspades    5.12   0   
##  6 metavelvet21  4.83   0   
##  7 metavelvet63  4.69   0   
##  8 metavelvet99 12.1    0   
##  9 ray21         6.20   0   
## 10 ray63         5.31   0   
## 11 ray99         5.31   0   
## 12 spades       18.4    3.35
## 13 trinity       7.00   0   
## 14 velvet21      7.58   0   
## 15 velvet63      4.79   0   
## 16 velvet99     11.9    0
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-8.png)<!-- -->

```r
make_boxplot("Total aligned length")
```

```
## # A tibble: 16 x 3
##    X2              mean median
##    <fct>          <dbl>  <dbl>
##  1 abyss21       23687. 28510 
##  2 abyss63       24264. 29344 
##  3 abyss99       22307. 27867 
##  4 megahit       27659. 29832 
##  5 metaspades    28961. 29829 
##  6 metavelvet21   5363.  2538 
##  7 metavelvet63   8408.  5183 
##  8 metavelvet99  11575.  8101 
##  9 ray21         29317. 29853 
## 10 ray63         28606. 29835 
## 11 ray99         28601. 29835 
## 12 spades       131003. 29873 
## 13 trinity       36969. 31436.
## 14 velvet21       5373.  2538 
## 15 velvet63       8485.  5217 
## 16 velvet99      11694.  7915
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-9.png)<!-- -->

```r
make_boxplot("misassemblies")
```

```
## # A tibble: 16 x 3
##    X2              mean median
##    <fct>          <dbl>  <dbl>
##  1 abyss21      0.00912      0
##  2 abyss63      0.00898      0
##  3 abyss99      0.0966       0
##  4 megahit      0.186        0
##  5 metaspades   0.349        0
##  6 metavelvet21 0.01         0
##  7 metavelvet63 0            0
##  8 metavelvet99 0.155        0
##  9 ray21        0.621        0
## 10 ray63        0.477        0
## 11 ray99        0.469        0
## 12 spades       4.73         0
## 13 trinity      1.13         0
## 14 velvet21     0.01         0
## 15 velvet63     0            0
## 16 velvet99     0.166        0
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-10.png)<!-- -->

```r
#make_boxplot("unaligned contigs")
make_boxplot("N50")
```

```
## # A tibble: 16 x 3
##    X2             mean median
##    <fct>         <dbl>  <dbl>
##  1 abyss21      11009.  5752 
##  2 abyss63       9622.  3385 
##  3 abyss99       7740.  1311 
##  4 megahit      16806. 22321 
##  5 metaspades   16498. 21586 
##  6 metavelvet21   649.   602 
##  7 metavelvet63   786.   682 
##  8 metavelvet99  1253.   732 
##  9 ray21        14863.  9286.
## 10 ray63        14523.  8450 
## 11 ray99        14525.  8450 
## 12 spades       14672. 10168.
## 13 trinity      13117.  7790 
## 14 velvet21       655.   598.
## 15 velvet63       792.   685 
## 16 velvet99      1235.   724
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-11.png)<!-- -->

```r
make_boxplot("N75")
```

```
## # A tibble: 16 x 3
##    X2             mean median
##    <fct>         <dbl>  <dbl>
##  1 abyss21       9331.  3087 
##  2 abyss63       7891.  1905 
##  3 abyss99       6749.   885 
##  4 megahit      13752.  4616.
##  5 metaspades   13081.  4740 
##  6 metavelvet21   575.   550 
##  7 metavelvet63   634.   574.
##  8 metavelvet99   930.   599 
##  9 ray21        11337.  3298.
## 10 ray63        11517.  2839 
## 11 ray99        11519.  2839 
## 12 spades       10810.  4634 
## 13 trinity       8391.  1706.
## 14 velvet21       576.   549 
## 15 velvet63       638.   574 
## 16 velvet99       933.   598
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-12.png)<!-- -->

```r
make_boxplot("L50")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21      32.5       2
##  2 abyss63      16.0       3
##  3 abyss99       6.32      3
##  4 megahit      47.2       1
##  5 metaspades   54.1       1
##  6 metavelvet21 32.4       4
##  7 metavelvet63 23.1       5
##  8 metavelvet99  9.37      5
##  9 ray21        20.8       2
## 10 ray63        17.9       2
## 11 ray99        17.9       2
## 12 spades       46.4       2
## 13 trinity      60.1       2
## 14 velvet21     32.6       4
## 15 velvet63     23.3       5
## 16 velvet99      9.38      5
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-13.png)<!-- -->

```r
make_boxplot("L75")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21       57.9      3
##  2 abyss63       29.6      5
##  3 abyss99       11.4      6
##  4 megahit       93.2      2
##  5 metaspades   105.       2
##  6 metavelvet21  53.4      6
##  7 metavelvet63  40.6      8
##  8 metavelvet99  16.1      8
##  9 ray21         42.3      3
## 10 ray63         36.2      3
## 11 ray99         36.2      3
## 12 spades       100.       3
## 13 trinity      118.       5
## 14 velvet21      53.6      6
## 15 velvet63      41.0      8
## 16 velvet99      16.1      8
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-14.png)<!-- -->

```r
make_boxplot("NA50")
```

```
## # A tibble: 16 x 3
##    X2             mean median
##    <fct>         <dbl>  <dbl>
##  1 abyss21      15552. 19288 
##  2 abyss63      13306.  7466 
##  3 abyss99       9868.  2256 
##  4 megahit      24922. 29821 
##  5 metaspades   24097. 28289 
##  6 metavelvet21   751.   666 
##  7 metavelvet63   881.   687 
##  8 metavelvet99  1420.   706 
##  9 ray21        18456. 22306 
## 10 ray63        18583. 22313 
## 11 ray99        18587. 22313 
## 12 spades       19237. 21999 
## 13 trinity      18204. 22300.
## 14 velvet21       778.   675 
## 15 velvet63       878.   703 
## 16 velvet99      1390.   686
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-15.png)<!-- -->

```r
make_boxplot("NA75")
```

```
## # A tibble: 16 x 3
##    X2             mean median
##    <fct>         <dbl>  <dbl>
##  1 abyss21      13513. 10223 
##  2 abyss63      11385.  4721 
##  3 abyss99       8800.  1277 
##  4 megahit      23024. 29821 
##  5 metaspades   21348. 27994 
##  6 metavelvet21   622.   576.
##  7 metavelvet63   684.   579 
##  8 metavelvet99  1045.   592.
##  9 ray21        15285. 10416.
## 10 ray63        16375. 13600 
## 11 ray99        16377. 13600 
## 12 spades       15877. 11043 
## 13 trinity      13684.  8573 
## 14 velvet21       627.   576.
## 15 velvet63       694.   582 
## 16 velvet99      1049.   598.
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-16.png)<!-- -->

```r
make_boxplot("LA50")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21       2.05      1
##  2 abyss63       2.65      2
##  3 abyss99       4.12      3
##  4 megahit       1.33      1
##  5 metaspades    1.44      1
##  6 metavelvet21  3.59      3
##  7 metavelvet63  4.27      4
##  8 metavelvet99  5.19      5
##  9 ray21         2.12      1
## 10 ray63         2.11      1
## 11 ray99         2.10      1
## 12 spades        5.64      1
## 13 trinity       3.02      1
## 14 velvet21      3.57      3
## 15 velvet63      4.41      4
## 16 velvet99      5.27      5
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-17.png)<!-- -->

```r
make_boxplot("LA75")
```

```
## # A tibble: 16 x 3
##    X2            mean median
##    <fct>        <dbl>  <dbl>
##  1 abyss21       3.12    2  
##  2 abyss63       3.98    3  
##  3 abyss99       7.07    6  
##  4 megahit       1.81    1  
##  5 metaspades    1.90    1  
##  6 metavelvet21  5.71    4.5
##  7 metavelvet63  6.43    6.5
##  8 metavelvet99  8.31    8  
##  9 ray21         3.27    2  
## 10 ray63         3.14    2  
## 11 ray99         3.15    2  
## 12 spades        7.11    2  
## 13 trinity       4.04    2  
## 14 velvet21      5.69    4  
## 15 velvet63      6.46    6  
## 16 velvet99      8.43    8
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-18.png)<!-- -->

```r
# genomic feature
x4 = xm %>% filter(Assembly == "genomic features") 
x4$value.1 = gsub("part", "", x4$value.1)
x5 = data.frame(str_split_fixed(x4$value.1,"\\+", 2), x4) 
x6 = data.frame(assembly = x5$X2, match = as.numeric(as.character(x5$X1.1)), mismatch = as.numeric(as.character(x5$X2.1)))

# need to add assay_type
  ggplot(x6, aes(assembly, 100*(match/49))) +
    geom_boxplot(aes(color = assembly)) +
    geom_jitter(alpha = .5, width = .1, aes(color = assembly)) +
    xlab("") +
    ylab("% features mapped") +
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = .5),
        axis.text = element_text(color = "black", angle = 90, hjust = 1)) 
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-2-19.png)<!-- -->

```r
ggsave("percent features mapped.pdf", width = 15, height = 10, units = "cm")
```

### 90% of the genome is covered by single contig


```r
library(tidyverse)
setwd("/Volumes/rony/drive/asm/covid19-Assembly/plots/")

#90%=29903*.9
contig = x3 %>% filter(Assembly == "Largest contig") %>% filter(value >=26912.7)
ggplot(contig, aes(X2)) +
    geom_histogram(stat = "count") +
    xlab("") +
    ylab("Count") +
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = .5),
        axis.text = element_text(color = "black", angle = 90, hjust = 1)) 
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

```r
ggsave("90_percent_genome_single_contig.pdf", width = 15, height = 10, units = "cm")
```

### correlation read vs genome


```r
library(tidyverse)
setwd("/Volumes/rony/drive/asm/covid19-Assembly/plots/")
r = read.table("~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/read_QC_matrix.txt")[,c(1,4)]
r = read.table("/Volumes/rony/drive/asm/covid19-Assembly/files/read_QC_matrix.txt")[,c(1,4)]

r2 = data.frame(id = str_split_fixed(r$V1, "_", 2), read = r$V4)
head(r2$id[,1])
```

```
## NULL
```

```r
#sum every two rows of PE data
r3 = data.frame(id = unique(r2$id.1), read = (rowsum(r2[,3], as.integer(gl(nrow(r2), 2, nrow(r2))))))
rx = xm %>% filter(Assembly == "Genome fraction (%)")
rx2 = inner_join(r3, rx, by = c("id" = "X1"))
#rxm = left_join(rx2, meta2, by = c("id" = "Run"))
rx2 %>% subset(is.na(rx2)) %>% dim()
```

```
## [1] 1115    9
```

```r
rx2 %>% subset(!is.na(rx2)) %>% dim()
```

```
## [1] 53128     9
```

```r
ggplot(rx2, aes(value, read/1e6)) +
    geom_point(aes(color = X2)) +
    geom_smooth(method='lm', formula= y~x) +
    xlab("Genome fraction (%)") +
    ylab("Number of reads (million)") + 
    facet_grid(Assay_Type~., scales = "free") +
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = .5),
        strip.background =element_rect(fill="white"),
        axis.text = element_text(color = "black", angle = 0, hjust = 1)) 
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
ggsave("readVsGenome.pdf", width = 12, height = 30, units = "cm")

rx3 = na.omit(rx2)
cor(rx3$value, rx3$read, method = "spearman")
```

```
## [1] 0.09714612
```

```r
#read dist; not needed
```


### sample to assembler


```r
library(tidyverse)
setwd("/Volumes/rony/drive/asm/covid19-Assembly/plots/")

rx2 %>% filter(Assembly == "Genome fraction (%)") %>% na.omit() %>%
    ggplot(aes(fct_reorder(id, read), X2, color = X2)) +
    geom_point(aes(size = value, alpha = .5)) +
    ylab("") +
    xlab("Samples are sorted by number of reads") + 
    theme(panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = .5),
        axis.text = element_text(color = "black", angle = 0, hjust = 1),
        axis.text.x=element_blank()) 
```

![](Assembly_analysis_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
ggsave("sampleVsassembler_color_bk.pdf", width = 30, height = 20, units = "cm")

# select 4 samples with high and low mean frac to plot in mvista
#high = SRR11578289 (97.4), SRR11597206 (94.4)
#low = SRR11828432, SRR11828424
#ref: MN908947.3

mvista = x3 %>% filter(Assembly == "Genome fraction (%)") %>% na.omit() %>%
    group_by(X1) %>%
    summarise(mean = mean(value)) %>% arrange(mean)
```


### align fasta to reference


```r
#https://www.biostars.org/p/110213/
# Build reference genome database
cd /projects/epigenomics3/temp/rislam/assembly/rajan/asm_pe/output_bioRxiv_100samples/ref/

gmap_build -D dir/ -d refgenome MN908947.3.fasta
# Alignment
gmap -D dir/refgenome/ -d refgenome -f samse -t 8 ../*fasta | samtools view -Shb - | samtools sort - alignment

samtools index aligment.bam
```


### dataset 


```r
library(tidyverse)
library(knitr)

c = read_tsv("~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/SraRunTable_COVID19_14.06.20.txt")

c2 = tibble(c$Platform, c$Run, c$SRA_Sample, c$Instrument, c$LibraryLayout, c$Assay_Type, c$LibrarySelection, c$LibrarySource, c$Organism, c$geo_loc_name, c$host, c$host_disease, c$Consent)

names(c2) <- gsub("c\\$", "", names(c2))
colSums(!is.na(c2))
```

```
##         Platform              Run       SRA_Sample       Instrument 
##            15007            15007            15007            15007 
##    LibraryLayout       Assay_Type LibrarySelection    LibrarySource 
##            15007            15007            15007            15007 
##         Organism     geo_loc_name             host     host_disease 
##            15007             5508             5373             5297 
##          Consent 
##            15007
```

```r
#
#write_csv(c2, "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/COVID19_14.06.20_metadata_final.csv")

# summarise metadata
#colnames(c2)

#good code example: https://uc-r.github.io/descriptives_categorical
table3 <- table(c2$Instrument, c2$Assay_Type)
table3 <- table( c2$Assay_Type, c2$LibrarySource)
table3 <- table(c2$Assay_Type, c2$LibrarySource, c2$LibraryLayout)
ftable(table3)
```

```
##                                      PAIRED SINGLE
##                                                   
## AMPLICON         GENOMIC                  0      0
##                  METAGENOMIC              0     41
##                  METATRANSCRIPTOMIC       0      0
##                  SYNTHETIC                0    136
##                  TRANSCRIPTOMIC           0      0
##                  VIRAL RNA             5254   6507
## OTHER            GENOMIC                  0      0
##                  METAGENOMIC              1      0
##                  METATRANSCRIPTOMIC       0      0
##                  SYNTHETIC                0      0
##                  TRANSCRIPTOMIC           0      0
##                  VIRAL RNA               68      0
## RNA-Seq          GENOMIC                  1      0
##                  METAGENOMIC              9      9
##                  METATRANSCRIPTOMIC       9      0
##                  SYNTHETIC                0      0
##                  TRANSCRIPTOMIC          16      3
##                  VIRAL RNA              682    654
## Targeted-Capture GENOMIC                  0      0
##                  METAGENOMIC              0      0
##                  METATRANSCRIPTOMIC       0      0
##                  SYNTHETIC                0      0
##                  TRANSCRIPTOMIC           0      0
##                  VIRAL RNA              194    241
## WGA              GENOMIC                  6      0
##                  METAGENOMIC              5      0
##                  METATRANSCRIPTOMIC       0      0
##                  SYNTHETIC                0      0
##                  TRANSCRIPTOMIC           0      0
##                  VIRAL RNA             1061      0
## WGS              GENOMIC                 11      0
##                  METAGENOMIC              6      0
##                  METATRANSCRIPTOMIC       0      0
##                  SYNTHETIC                0      0
##                  TRANSCRIPTOMIC           0      0
##                  VIRAL RNA               93      0
```

```r
# will add table paper
df = c2 %>% 
    group_by(Assay_Type, LibrarySource, LibraryLayout) %>%
    tally() 

kable(df, caption = "Summary of all data. This table will add table paper")
```



Table: Summary of all data. This table will add table paper

|Assay_Type       |LibrarySource      |LibraryLayout |    n|
|:----------------|:------------------|:-------------|----:|
|AMPLICON         |METAGENOMIC        |SINGLE        |   41|
|AMPLICON         |SYNTHETIC          |SINGLE        |  136|
|AMPLICON         |VIRAL RNA          |PAIRED        | 5254|
|AMPLICON         |VIRAL RNA          |SINGLE        | 6507|
|OTHER            |METAGENOMIC        |PAIRED        |    1|
|OTHER            |VIRAL RNA          |PAIRED        |   68|
|RNA-Seq          |GENOMIC            |PAIRED        |    1|
|RNA-Seq          |METAGENOMIC        |PAIRED        |    9|
|RNA-Seq          |METAGENOMIC        |SINGLE        |    9|
|RNA-Seq          |METATRANSCRIPTOMIC |PAIRED        |    9|
|RNA-Seq          |TRANSCRIPTOMIC     |PAIRED        |   16|
|RNA-Seq          |TRANSCRIPTOMIC     |SINGLE        |    3|
|RNA-Seq          |VIRAL RNA          |PAIRED        |  682|
|RNA-Seq          |VIRAL RNA          |SINGLE        |  654|
|Targeted-Capture |VIRAL RNA          |PAIRED        |  194|
|Targeted-Capture |VIRAL RNA          |SINGLE        |  241|
|WGA              |GENOMIC            |PAIRED        |    6|
|WGA              |METAGENOMIC        |PAIRED        |    5|
|WGA              |VIRAL RNA          |PAIRED        | 1061|
|WGS              |GENOMIC            |PAIRED        |   11|
|WGS              |METAGENOMIC        |PAIRED        |    6|
|WGS              |VIRAL RNA          |PAIRED        |   93|

```r
#write_tsv(df, "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/summary_data.tsv")

df2 = df %>% filter(LibrarySource == "VIRAL RNA")
kable(df2, caption = "Summary of VIRAL RNA data")
```



Table: Summary of VIRAL RNA data

|Assay_Type       |LibrarySource |LibraryLayout |    n|
|:----------------|:-------------|:-------------|----:|
|AMPLICON         |VIRAL RNA     |PAIRED        | 5254|
|AMPLICON         |VIRAL RNA     |SINGLE        | 6507|
|OTHER            |VIRAL RNA     |PAIRED        |   68|
|RNA-Seq          |VIRAL RNA     |PAIRED        |  682|
|RNA-Seq          |VIRAL RNA     |SINGLE        |  654|
|Targeted-Capture |VIRAL RNA     |PAIRED        |  194|
|Targeted-Capture |VIRAL RNA     |SINGLE        |  241|
|WGA              |VIRAL RNA     |PAIRED        | 1061|
|WGS              |VIRAL RNA     |PAIRED        |   93|

```r
##subsample main paper
#PE
set.seed(2020)
a1 = c2 %>% filter(LibraryLayout == "PAIRED" & Assay_Type == "AMPLICON" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "PE: AMPLICON of VIRAL RNA") %>% sample_n(100)

set.seed(2020)
a2 = c2 %>% filter(LibraryLayout == "PAIRED" & Assay_Type == "OTHER" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "PE: OTHER of VIRAL RNA")

set.seed(2020)
a3 = c2 %>% filter(LibraryLayout == "PAIRED" & Assay_Type == "RNA-Seq" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "PE: RNA-Seq of VIRAL RNA") %>% sample_n(100)

set.seed(2020)
a4 = c2 %>% filter(LibraryLayout == "PAIRED" & Assay_Type == "Targeted-Capture" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "PE: Targeted-Capture of VIRAL RNA") %>% sample_n(100)

set.seed(2020)
a5 = c2 %>% filter(LibraryLayout == "PAIRED" & Assay_Type == "WGA" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "PE: WGA of VIRAL RNA") %>% sample_n(100)

set.seed(2020)
a6 = c2 %>% filter(LibraryLayout == "PAIRED" & Assay_Type == "WGS" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "PE: WGS of VIRAL RNA")

#SE
set.seed(2020)
b1 = c2 %>% filter(LibraryLayout == "SINGLE" & Assay_Type == "AMPLICON" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "SE: AMPLICON of VIRAL RNA") %>% sample_n(100)

set.seed(2020)
b2 = c2 %>% filter(LibraryLayout == "SINGLE" & Assay_Type == "RNA-Seq" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "SE: RNA-Seq of VIRAL RNA") %>% sample_n(100)

set.seed(2020)
b3 = c2 %>% filter(LibraryLayout == "SINGLE" & Assay_Type == "Targeted-Capture" & LibrarySource == "VIRAL RNA") %>%
    mutate(LibType = "SE: Targeted-Capture of VIRAL RNA") %>% sample_n(100)

#dataset
ab_pe = rbind(a1, a2, a3, a4, a5, a6)
ab_se = rbind(b1, b2, b3)
#write_csv(ab_pe, "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/PE_561samples_final.csv")
#write_csv(ab_se, "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/SE_300samples_final.csv")
#write.table(ab_pe$Run, "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/PE_561samples_final_561runs.txt", col.names = F, row.names = F, quote = F)
#write.table(ab_se$Run, "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/SE_300samples_final_300runs.txt", col.names = F, row.names = F, quote = F)
#write_csv(a1, "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/PE_100samples_amplicon_bioRxiv.csv")
#write.table(a1[,2], "~/Gdrive_tutorial_edits/Assembly_COVID19/covid19-Assembly/files/PE_100samples_amplicon_bioRxiv_100runs.txt", col.names = F, row.names = F, quote = F)

# test ram and cpu
set.seed(2020)
#20 samples will be selected from amplicon with similar depth
#ram = a1 %>% sample_n(20)

ab = rbind(ab_pe, ab_se)
ab2 = ab %>% group_by(LibType) %>% tally()
kable(ab2, caption = "List of total 9 different categories. Maximum 100 samples are randomly selected")
```



Table: List of total 9 different categories. Maximum 100 samples are randomly selected

|LibType                           |   n|
|:---------------------------------|---:|
|PE: AMPLICON of VIRAL RNA         | 100|
|PE: OTHER of VIRAL RNA            |  68|
|PE: RNA-Seq of VIRAL RNA          | 100|
|PE: Targeted-Capture of VIRAL RNA | 100|
|PE: WGA of VIRAL RNA              | 100|
|PE: WGS of VIRAL RNA              |  93|
|SE: AMPLICON of VIRAL RNA         | 100|
|SE: RNA-Seq of VIRAL RNA          | 100|
|SE: Targeted-Capture of VIRAL RNA | 100|

