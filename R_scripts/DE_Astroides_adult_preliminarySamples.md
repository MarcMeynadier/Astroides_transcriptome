DE_Astroides_adult_preliminarySamples
================
Marc Meynadier
6/3/2022

``` r
# Differential expression on Kallisto data 

# Preliminary samples - 2016 dataset

# Packages and dependence
packageCheckClassic <- function(x){
  # 
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}

packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','pheatmap','RColorBrewer','genefilter','gplots','vegan','dplyr'))
```

    ## Le chargement a nécessité le package : DESeq2

    ## Le chargement a nécessité le package : S4Vectors

    ## Warning: le package 'S4Vectors' a été compilé avec la version R 4.1.3

    ## Le chargement a nécessité le package : stats4

    ## Le chargement a nécessité le package : BiocGenerics

    ## 
    ## Attachement du package : 'BiocGenerics'

    ## Les objets suivants sont masqués depuis 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## Les objets suivants sont masqués depuis 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## 
    ## Attachement du package : 'S4Vectors'

    ## Les objets suivants sont masqués depuis 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Le chargement a nécessité le package : IRanges

    ## Le chargement a nécessité le package : GenomicRanges

    ## Le chargement a nécessité le package : GenomeInfoDb

    ## Le chargement a nécessité le package : SummarizedExperiment

    ## Le chargement a nécessité le package : MatrixGenerics

    ## Le chargement a nécessité le package : matrixStats

    ## 
    ## Attachement du package : 'MatrixGenerics'

    ## Les objets suivants sont masqués depuis 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Le chargement a nécessité le package : Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attachement du package : 'Biobase'

    ## L'objet suivant est masqué depuis 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## Les objets suivants sont masqués depuis 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## Le chargement a nécessité le package : devtools

    ## Le chargement a nécessité le package : usethis

    ## Le chargement a nécessité le package : BiocManager

    ## Bioconductor version '3.14' is out-of-date; the current release version '3.15'
    ##   is available with R version '4.2'; see https://bioconductor.org/install

    ## 
    ## Attachement du package : 'BiocManager'

    ## L'objet suivant est masqué depuis 'package:devtools':
    ## 
    ##     install

    ## Le chargement a nécessité le package : ggplot2

    ## Le chargement a nécessité le package : ggrepel

    ## Le chargement a nécessité le package : markdown

    ## Le chargement a nécessité le package : pheatmap

    ## Le chargement a nécessité le package : RColorBrewer

    ## Le chargement a nécessité le package : genefilter

    ## 
    ## Attachement du package : 'genefilter'

    ## Les objets suivants sont masqués depuis 'package:MatrixGenerics':
    ## 
    ##     rowSds, rowVars

    ## Les objets suivants sont masqués depuis 'package:matrixStats':
    ## 
    ##     rowSds, rowVars

    ## Le chargement a nécessité le package : gplots

    ## 
    ## Attachement du package : 'gplots'

    ## L'objet suivant est masqué depuis 'package:IRanges':
    ## 
    ##     space

    ## L'objet suivant est masqué depuis 'package:S4Vectors':
    ## 
    ##     space

    ## L'objet suivant est masqué depuis 'package:stats':
    ## 
    ##     lowess

    ## Le chargement a nécessité le package : vegan

    ## Le chargement a nécessité le package : permute

    ## 
    ## Attachement du package : 'permute'

    ## L'objet suivant est masqué depuis 'package:devtools':
    ## 
    ##     check

    ## Le chargement a nécessité le package : lattice

    ## This is vegan 2.6-2

    ## Le chargement a nécessité le package : dplyr

    ## 
    ## Attachement du package : 'dplyr'

    ## L'objet suivant est masqué depuis 'package:Biobase':
    ## 
    ##     combine

    ## L'objet suivant est masqué depuis 'package:matrixStats':
    ## 
    ##     count

    ## Les objets suivants sont masqués depuis 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## L'objet suivant est masqué depuis 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## Les objets suivants sont masqués depuis 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## Les objets suivants sont masqués depuis 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## Les objets suivants sont masqués depuis 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## Les objets suivants sont masqués depuis 'package:stats':
    ## 
    ##     filter, lag

    ## Les objets suivants sont masqués depuis 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
```

    ## Skipping install of 'ggvenn' from a github remote, the SHA1 (b7ff54ba) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library('ggvenn')
```

    ## Le chargement a nécessité le package : grid

``` r
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
```

    ## Registered S3 methods overwritten by 'ggalt':
    ##   method                  from   
    ##   grid.draw.absoluteGrob  ggplot2
    ##   grobHeight.absoluteGrob ggplot2
    ##   grobWidth.absoluteGrob  ggplot2
    ##   grobX.absoluteGrob      ggplot2
    ##   grobY.absoluteGrob      ggplot2

``` r
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
```

    ## ℹ SHA-1 hash of file is 015fc0457e61e3e93a903e69a24d96d2dac7b9fb

``` r
# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
samples<-read.table('tximport_design_preliminarySamples.txt',header=T)
tx2gene<-read.table('tx2gene_adultTranscriptome',header=T)
candidateGenes<-read.csv('candidateGenes.csv',header=T,sep=',')
scriptPath <- sub("/[^/]+$", "", scriptPath)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/6_kallisto/adultTranscriptome/adult/1_preliminarySamples'
outputPath<-paste(scriptPath,'/output/DESeq2/adultTranscriptome/adult/1_preliminarySamples/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
files<-paste0(samples$sample,'.tsv')
names(files)<-samples$sample
txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
```

    ## Note: importing `abundance.h5` is typically faster than `abundance.tsv`

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 
    ## transcripts missing from tx2gene: 1
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
names(txi)
```

    ## [1] "abundance"           "counts"              "length"             
    ## [4] "countsFromAbundance"

``` r
head(txi$counts)
```

    ##                    abundance_adult_nov2016_gm_gm_pre_14283X1_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                   828.000
    ## TRINITY_DN0_c0_g2                                                   237.000
    ## TRINITY_DN0_c1_g1                                                    61.000
    ## TRINITY_DN1_c0_g1                                                  3550.000
    ## TRINITY_DN1_c1_g1                                                     0.000
    ## TRINITY_DN10_c0_g1                                                  393.128
    ##                    abundance_adult_nov2016_gm_gm_pre_14283X2_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                   857.000
    ## TRINITY_DN0_c0_g2                                                   243.119
    ## TRINITY_DN0_c1_g1                                                    74.000
    ## TRINITY_DN1_c0_g1                                                  3361.000
    ## TRINITY_DN1_c1_g1                                                     0.000
    ## TRINITY_DN10_c0_g1                                                  387.428
    ##                    abundance_adult_nov2016_gm_gm_pre_14283X3_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                  675.0000
    ## TRINITY_DN0_c0_g2                                                  163.1960
    ## TRINITY_DN0_c1_g1                                                   44.2649
    ## TRINITY_DN1_c0_g1                                                 6772.0000
    ## TRINITY_DN1_c1_g1                                                    0.0000
    ## TRINITY_DN10_c0_g1                                                 470.3850
    ##                    abundance_adult_nov2016_gm_gm_pre_14283X4_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                   977.000
    ## TRINITY_DN0_c0_g2                                                   229.763
    ## TRINITY_DN0_c1_g1                                                    51.000
    ## TRINITY_DN1_c0_g1                                                  4706.000
    ## TRINITY_DN1_c1_g1                                                     0.000
    ## TRINITY_DN10_c0_g1                                                  224.756
    ##                    abundance_adult_nov2016_pv_pv_pre_14283X10_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                   727.0000
    ## TRINITY_DN0_c0_g2                                                   171.2530
    ## TRINITY_DN0_c1_g1                                                    51.2069
    ## TRINITY_DN1_c0_g1                                                  3797.0000
    ## TRINITY_DN1_c1_g1                                                     1.0000
    ## TRINITY_DN10_c0_g1                                                  518.4680
    ##                    abundance_adult_nov2016_pv_pv_pre_14283X8_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                   965.000
    ## TRINITY_DN0_c0_g2                                                   214.888
    ## TRINITY_DN0_c1_g1                                                    44.000
    ## TRINITY_DN1_c0_g1                                                  4569.000
    ## TRINITY_DN1_c1_g1                                                     0.000
    ## TRINITY_DN10_c0_g1                                                  407.287
    ##                    abundance_adult_nov2016_pv_pv_pre_14283X9_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                  522.0000
    ## TRINITY_DN0_c0_g2                                                  156.8880
    ## TRINITY_DN0_c1_g1                                                   56.3134
    ## TRINITY_DN1_c0_g1                                                 4562.0000
    ## TRINITY_DN1_c1_g1                                                    1.0000
    ## TRINITY_DN10_c0_g1                                                 585.4630
    ##                    abundance_adult_nov2016_sa_sa_pre_14283X5_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                  1169.000
    ## TRINITY_DN0_c0_g2                                                   262.969
    ## TRINITY_DN0_c1_g1                                                    26.000
    ## TRINITY_DN1_c0_g1                                                  3466.000
    ## TRINITY_DN1_c1_g1                                                     0.000
    ## TRINITY_DN10_c0_g1                                                  490.226
    ##                    abundance_adult_nov2016_sa_sa_pre_14283X6_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                   617.000
    ## TRINITY_DN0_c0_g2                                                   404.000
    ## TRINITY_DN0_c1_g1                                                     6.000
    ## TRINITY_DN1_c0_g1                                                  5516.000
    ## TRINITY_DN1_c1_g1                                                     0.000
    ## TRINITY_DN10_c0_g1                                                  340.611
    ##                    abundance_adult_nov2016_sa_sa_pre_14283X7_paired_trimmed
    ## TRINITY_DN0_c0_g1                                                   800.000
    ## TRINITY_DN0_c0_g2                                                   236.890
    ## TRINITY_DN0_c1_g1                                                    63.000
    ## TRINITY_DN1_c0_g1                                                  4154.000
    ## TRINITY_DN1_c1_g1                                                     0.000
    ## TRINITY_DN10_c0_g1                                                  554.357

``` r
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

``` r
# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

# Differential expression analysis
dds<-DESeq(dds)
```

    ## estimating size factors
    ## using 'avgTxLength' from assays(dds), correcting for library size
    ## estimating dispersions
    ## gene-wise dispersion estimates
    ## mean-dispersion relationship
    ## final dispersion estimates
    ## fitting model and testing

``` r
cbind(resultsNames(dds))
```

    ##      [,1]           
    ## [1,] "Intercept"    
    ## [2,] "site_pv_vs_gm"
    ## [3,] "site_sa_vs_gm"

``` r
res_pv_gm<-results(dds, contrast=c("site","pv","gm"), alpha = 0.05)
res_sa_gm<-results(dds, contrast=c("site","sa","gm"), alpha = 0.05)
res_pv_sa<-results(dds, contrast=c("site","pv","sa"), alpha = 0.05)

# Exploring the results

# Results pv VS gm

#MA-plot
resLFC = lfcShrink(dds, contrast=c("site","pv","gm"), 
                   type="ashr")
```

    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041

``` r
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes\nPreliminary samples : pv VS gm")
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
EnhancedVolcano(data.frame(res_pv_gm), lab = rownames(data.frame(res_pv_gm)), x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "Preliminary samples : pv VS gm",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_gm), ' variables'),
                    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
# Results sa VS gm

#MA-plot
resLFC = lfcShrink(dds, contrast=c("site","sa","gm"), 
                   type="ashr")
```

    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041

``` r
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes\nPreliminary samples : sa VS gm")
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
# Volcano plot
EnhancedVolcano(data.frame(res_sa_gm), lab = rownames(data.frame(res_sa_gm)), x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "Preliminary samples : sa VS gm",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_sa_gm), ' variables'),
                    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->

``` r
# Results pv VS sa

#MA-plot
resLFC = lfcShrink(dds, contrast=c("site","pv","sa"), 
                   type="ashr")
```

    ## using 'ashr' for LFC shrinkage. If used in published research, please cite:
    ##     Stephens, M. (2016) False discovery rates: a new deal. Biostatistics, 18:2.
    ##     https://doi.org/10.1093/biostatistics/kxw041

``` r
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes\nPreliminary samples : pv VS sa")
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->

``` r
# Volcano plot
EnhancedVolcano(data.frame(res_pv_sa), lab = rownames(data.frame(res_pv_sa)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Preliminary samples : pv VS sa",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_sa), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->

``` r
# Principal Component Analysis

vsd = vst(dds,blind=T)

pcaData = plotPCA(vsd, intgroup="site", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of adult corals", subtitle = "nov2016 dataset") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-7.png)<!-- -->

``` r
# Venn diagramm 
resOrdered_pv_gm <- res_pv_gm[order(res_pv_gm$padj),]
resOrderedDF_pv_gm <- as.data.frame(resOrdered_pv_gm)
resOrderedDF_pv_gm_venn <- filter(resOrderedDF_pv_gm,padj < 0.05)
resOrderedDF_pv_gm_venn <- list(rownames(resOrderedDF_pv_gm_venn))
resOrderedDF_pv_gm_venn <- unlist(resOrderedDF_pv_gm_venn)

resOrdered_sa_gm <- res_sa_gm[order(res_sa_gm$padj),]
resOrderedDF_sa_gm <- as.data.frame(resOrdered_sa_gm)
resOrderedDF_sa_gm_venn <- filter(resOrderedDF_sa_gm,padj < 0.05)
resOrderedDF_sa_gm_venn <- list(rownames(resOrderedDF_sa_gm_venn))
resOrderedDF_sa_gm_venn <- unlist(resOrderedDF_sa_gm_venn)

resOrdered_pv_sa <- res_pv_sa[order(res_pv_sa$padj),]
resOrderedDF_pv_sa <- as.data.frame(resOrdered_pv_sa)
resOrderedDF_pv_sa_venn <- filter(resOrderedDF_pv_sa,padj < 0.05)
resOrderedDF_pv_sa_venn <- list(rownames(resOrderedDF_pv_sa_venn))
resOrderedDF_pv_sa_venn <- unlist(resOrderedDF_pv_sa_venn)

x = list('pv VS sa' = resOrderedDF_pv_sa_venn, 'sa VS gm' = resOrderedDF_sa_gm_venn, 'pv VS gm' = resOrderedDF_pv_gm_venn)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-8.png)<!-- -->

``` r
# Candidate genes heatmap

listGenes <- candidateGenes$genes
listGenes2 <- which(rownames(vsd) %in% listGenes)
index <- which(listGenes %in% rownames(vsd))
candidateGenes2 <- candidateGenes[index, ] 
listProt <- candidateGenes2$pfam_annotation
listGenes3 <- candidateGenes2$genes

vsdCandidate <- vsd[listGenes3, ]

labColName <- c('gm','gm','gm','gm','pv','pv','pv','sa','sa','sa')
colnames(vsdCandidate) <- labColName
rownames(vsdCandidate) <- listProt

topVarGenesVsd <- head(order(rowVars(assay(vsdCandidate)), decreasing=TRUE), 50 )
heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",
          key.title = "",
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7,density.info="none",
          xlab="sampling sites",ylab="proteins associated to genes",Colv=NA,margins = c(4, 7))
```

    ## Warning in heatmap.2(assay(vsdCandidate)[topVarGenesVsd, ], trace = "none", :
    ## Discrepancy: Colv is FALSE, while dendrogram is `both'. Omitting column
    ## dendogram.

``` r
main='Differential expression of 50 most expressed candidates genes\n\nPreliminary samples'
title(main, cex.main = 0.7)
```

![](DE_Astroides_adult_preliminarySamples_files/figure-gfm/unnamed-chunk-1-9.png)<!-- -->

``` r
# Inferences statistics
count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site, method="euclidian")
```

    ## 'adonis' will be deprecated: use 'adonis2' instead

    ## $aov.tab
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## site       2     55701   27850   1.722 0.32975  0.001 ***
    ## Residuals  7    113216   16174         0.67025           
    ## Total      9    168917                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $call
    ## adonis(formula = dist_tab_assay ~ site, data = samples, method = "euclidian")
    ## 
    ## $coefficients
    ## NULL
    ## 
    ## $coef.sites
    ##                  [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
    ## (Intercept) 173.32456 184.02046 184.16643 173.75365 175.37706 170.59302
    ## site1       -36.06244 -42.71170 -38.58996 -37.26829  25.98829  23.77230
    ## site2        18.65456  20.40694  12.26533  13.63107 -59.03986 -49.41587
    ##                  [,7]      [,8]      [,9]     [,10]
    ## (Intercept) 163.60748 172.11040 174.67845 167.06240
    ## site1        25.82909  29.78205  29.25691  30.81448
    ## site2       -50.58663  22.20415  30.30689  26.66796
    ## 
    ## $f.perms
    ##              [,1]
    ##    [1,] 1.1930764
    ##    [2,] 1.0067114
    ##    [3,] 1.1236261
    ##    [4,] 0.8588821
    ##    [5,] 0.8357488
    ##    [6,] 1.1116440
    ##    [7,] 0.7921230
    ##    [8,] 0.9523803
    ##    [9,] 1.0509783
    ##   [10,] 1.1305179
    ##   [11,] 0.8066922
    ##   [12,] 1.0638217
    ##   [13,] 0.9859505
    ##   [14,] 1.1785664
    ##   [15,] 1.0440657
    ##   [16,] 1.2022336
    ##   [17,] 0.9898347
    ##   [18,] 0.8936216
    ##   [19,] 1.0449349
    ##   [20,] 1.0440657
    ##   [21,] 0.9817481
    ##   [22,] 1.1102552
    ##   [23,] 1.0381492
    ##   [24,] 0.8366617
    ##   [25,] 0.9867058
    ##   [26,] 0.9599648
    ##   [27,] 0.9725166
    ##   [28,] 0.9482863
    ##   [29,] 1.1414865
    ##   [30,] 1.1510094
    ##   [31,] 1.1004839
    ##   [32,] 1.1223893
    ##   [33,] 1.0305753
    ##   [34,] 0.9282987
    ##   [35,] 0.9134379
    ##   [36,] 1.1697229
    ##   [37,] 0.7783088
    ##   [38,] 1.0290787
    ##   [39,] 0.9351232
    ##   [40,] 0.7540891
    ##   [41,] 0.9369821
    ##   [42,] 0.8898014
    ##   [43,] 0.9313278
    ##   [44,] 1.1247513
    ##   [45,] 0.9510644
    ##   [46,] 0.8490154
    ##   [47,] 0.9827484
    ##   [48,] 1.0218927
    ##   [49,] 1.2149461
    ##   [50,] 0.8709922
    ##   [51,] 0.7951333
    ##   [52,] 1.0869203
    ##   [53,] 1.4788897
    ##   [54,] 0.8081432
    ##   [55,] 0.9725166
    ##   [56,] 1.0638574
    ##   [57,] 1.1089728
    ##   [58,] 1.0578503
    ##   [59,] 0.9118340
    ##   [60,] 1.0694514
    ##   [61,] 0.9306455
    ##   [62,] 1.1379190
    ##   [63,] 0.9099921
    ##   [64,] 1.0845869
    ##   [65,] 0.9221634
    ##   [66,] 1.0550483
    ##   [67,] 1.1290570
    ##   [68,] 0.8039192
    ##   [69,] 0.8378759
    ##   [70,] 1.0488246
    ##   [71,] 0.8216925
    ##   [72,] 1.0797307
    ##   [73,] 1.0857646
    ##   [74,] 0.9282987
    ##   [75,] 0.8697745
    ##   [76,] 0.9192675
    ##   [77,] 1.0102917
    ##   [78,] 0.9852926
    ##   [79,] 0.9406028
    ##   [80,] 0.9913084
    ##   [81,] 1.1607894
    ##   [82,] 0.8304142
    ##   [83,] 0.8840747
    ##   [84,] 0.9945006
    ##   [85,] 1.2654087
    ##   [86,] 1.1065117
    ##   [87,] 0.9941660
    ##   [88,] 1.0106091
    ##   [89,] 0.8547079
    ##   [90,] 0.7987765
    ##   [91,] 0.9008549
    ##   [92,] 0.8824936
    ##   [93,] 0.8617930
    ##   [94,] 1.0400463
    ##   [95,] 0.8757593
    ##   [96,] 0.8480838
    ##   [97,] 1.1248285
    ##   [98,] 1.3130098
    ##   [99,] 0.8774938
    ##  [100,] 1.3555407
    ##  [101,] 0.8286825
    ##  [102,] 0.8423008
    ##  [103,] 0.9884835
    ##  [104,] 1.0649420
    ##  [105,] 0.9928384
    ##  [106,] 1.1064433
    ##  [107,] 1.2389125
    ##  [108,] 0.7960225
    ##  [109,] 1.0014831
    ##  [110,] 1.0376123
    ##  [111,] 0.9764786
    ##  [112,] 0.8490154
    ##  [113,] 0.8365371
    ##  [114,] 1.1254892
    ##  [115,] 0.8665142
    ##  [116,] 0.9406028
    ##  [117,] 1.0188752
    ##  [118,] 1.1037622
    ##  [119,] 1.0063992
    ##  [120,] 0.9770904
    ##  [121,] 0.9031565
    ##  [122,] 0.8881486
    ##  [123,] 1.0317230
    ##  [124,] 0.8674496
    ##  [125,] 1.1561520
    ##  [126,] 0.8877536
    ##  [127,] 0.8840747
    ##  [128,] 0.8022045
    ##  [129,] 0.9304281
    ##  [130,] 1.0970970
    ##  [131,] 1.1247513
    ##  [132,] 0.7675518
    ##  [133,] 1.1288975
    ##  [134,] 0.9607405
    ##  [135,] 1.0798853
    ##  [136,] 1.0768972
    ##  [137,] 0.9217873
    ##  [138,] 1.1386574
    ##  [139,] 0.8918536
    ##  [140,] 0.9465700
    ##  [141,] 1.0211861
    ##  [142,] 0.9803131
    ##  [143,] 0.8582562
    ##  [144,] 0.9938178
    ##  [145,] 0.8530780
    ##  [146,] 1.0158603
    ##  [147,] 0.7398403
    ##  [148,] 1.1093801
    ##  [149,] 1.2794941
    ##  [150,] 1.1290570
    ##  [151,] 1.0218927
    ##  [152,] 0.8222756
    ##  [153,] 0.8945047
    ##  [154,] 0.9467812
    ##  [155,] 0.8378759
    ##  [156,] 1.0353567
    ##  [157,] 1.1266987
    ##  [158,] 0.8288004
    ##  [159,] 1.0806734
    ##  [160,] 1.1697229
    ##  [161,] 0.9094614
    ##  [162,] 1.0687103
    ##  [163,] 1.0178981
    ##  [164,] 0.8534052
    ##  [165,] 0.9491554
    ##  [166,] 1.0998046
    ##  [167,] 0.8385139
    ##  [168,] 0.9615020
    ##  [169,] 0.9580785
    ##  [170,] 1.0509783
    ##  [171,] 1.1067805
    ##  [172,] 1.0870736
    ##  [173,] 0.8295287
    ##  [174,] 0.9008634
    ##  [175,] 0.7782556
    ##  [176,] 0.9635467
    ##  [177,] 0.8209148
    ##  [178,] 0.8304142
    ##  [179,] 0.9318107
    ##  [180,] 1.0493231
    ##  [181,] 1.4451137
    ##  [182,] 1.0793743
    ##  [183,] 1.0376471
    ##  [184,] 1.0725983
    ##  [185,] 1.0764814
    ##  [186,] 0.8779192
    ##  [187,] 0.9134460
    ##  [188,] 1.0024654
    ##  [189,] 0.8931808
    ##  [190,] 1.0564406
    ##  [191,] 0.9530983
    ##  [192,] 0.9905325
    ##  [193,] 1.1330625
    ##  [194,] 0.9987543
    ##  [195,] 0.8677659
    ##  [196,] 1.1373092
    ##  [197,] 1.0261235
    ##  [198,] 1.0986444
    ##  [199,] 1.1396027
    ##  [200,] 0.8501382
    ##  [201,] 0.8969267
    ##  [202,] 1.0256205
    ##  [203,] 0.9923852
    ##  [204,] 0.9600742
    ##  [205,] 1.0057378
    ##  [206,] 0.8254301
    ##  [207,] 1.0620931
    ##  [208,] 0.8472885
    ##  [209,] 0.9878450
    ##  [210,] 0.9588629
    ##  [211,] 1.2860673
    ##  [212,] 0.8428152
    ##  [213,] 1.1199874
    ##  [214,] 0.9095234
    ##  [215,] 0.9104471
    ##  [216,] 1.0733738
    ##  [217,] 1.0191287
    ##  [218,] 1.2038557
    ##  [219,] 1.0739108
    ##  [220,] 1.0045405
    ##  [221,] 1.0775650
    ##  [222,] 0.9062112
    ##  [223,] 0.9094571
    ##  [224,] 1.0649924
    ##  [225,] 1.1223285
    ##  [226,] 1.0031994
    ##  [227,] 0.9661499
    ##  [228,] 1.0123041
    ##  [229,] 1.0405876
    ##  [230,] 1.0339656
    ##  [231,] 0.8230913
    ##  [232,] 1.0997155
    ##  [233,] 0.9635491
    ##  [234,] 1.4206564
    ##  [235,] 0.9965850
    ##  [236,] 1.1414865
    ##  [237,] 1.0886238
    ##  [238,] 1.1319114
    ##  [239,] 0.9362413
    ##  [240,] 0.8825076
    ##  [241,] 1.0426477
    ##  [242,] 0.9157333
    ##  [243,] 1.3331708
    ##  [244,] 1.0105664
    ##  [245,] 1.1865208
    ##  [246,] 0.9815266
    ##  [247,] 0.9788386
    ##  [248,] 0.8813392
    ##  [249,] 0.9658031
    ##  [250,] 1.0863082
    ##  [251,] 1.0024654
    ##  [252,] 0.9728224
    ##  [253,] 1.0358185
    ##  [254,] 0.9175692
    ##  [255,] 1.0282929
    ##  [256,] 0.8859677
    ##  [257,] 0.8313344
    ##  [258,] 1.0488246
    ##  [259,] 1.1370994
    ##  [260,] 1.1220130
    ##  [261,] 1.0400463
    ##  [262,] 0.7684228
    ##  [263,] 0.7800117
    ##  [264,] 1.1715330
    ##  [265,] 1.1247619
    ##  [266,] 0.8081909
    ##  [267,] 1.0385767
    ##  [268,] 0.9947389
    ##  [269,] 0.8800600
    ##  [270,] 1.0730549
    ##  [271,] 0.8859677
    ##  [272,] 1.1486772
    ##  [273,] 1.1643890
    ##  [274,] 0.9248986
    ##  [275,] 1.0141373
    ##  [276,] 0.9104471
    ##  [277,] 1.0280507
    ##  [278,] 1.0048754
    ##  [279,] 1.3792216
    ##  [280,] 0.8212487
    ##  [281,] 0.9925005
    ##  [282,] 1.0513241
    ##  [283,] 0.9740930
    ##  [284,] 0.9080173
    ##  [285,] 0.9465700
    ##  [286,] 0.9770443
    ##  [287,] 1.0511023
    ##  [288,] 0.9570790
    ##  [289,] 1.1467138
    ##  [290,] 1.1661902
    ##  [291,] 1.0381448
    ##  [292,] 1.0205491
    ##  [293,] 0.8483318
    ##  [294,] 0.8588821
    ##  [295,] 0.9799787
    ##  [296,] 1.1316020
    ##  [297,] 0.9901672
    ##  [298,] 1.0114224
    ##  [299,] 1.1229832
    ##  [300,] 0.9822325
    ##  [301,] 1.1278102
    ##  [302,] 1.0849909
    ##  [303,] 0.9938178
    ##  [304,] 0.9021073
    ##  [305,] 1.0611598
    ##  [306,] 1.0563017
    ##  [307,] 1.0160626
    ##  [308,] 0.9327114
    ##  [309,] 1.2677004
    ##  [310,] 0.9236968
    ##  [311,] 0.8386403
    ##  [312,] 0.8067463
    ##  [313,] 1.0080378
    ##  [314,] 1.2161462
    ##  [315,] 0.8887313
    ##  [316,] 1.1758543
    ##  [317,] 1.0268687
    ##  [318,] 0.9719667
    ##  [319,] 1.0708762
    ##  [320,] 1.0445303
    ##  [321,] 0.9040711
    ##  [322,] 0.9320734
    ##  [323,] 0.8917513
    ##  [324,] 0.9178154
    ##  [325,] 1.1979156
    ##  [326,] 0.7705444
    ##  [327,] 0.8775495
    ##  [328,] 0.9914357
    ##  [329,] 0.7741013
    ##  [330,] 1.1498174
    ##  [331,] 1.0531998
    ##  [332,] 0.9709668
    ##  [333,] 0.8295287
    ##  [334,] 0.9111871
    ##  [335,] 0.9682071
    ##  [336,] 0.9801847
    ##  [337,] 1.2040256
    ##  [338,] 0.9891640
    ##  [339,] 0.9963599
    ##  [340,] 1.0135141
    ##  [341,] 1.1275790
    ##  [342,] 0.9868017
    ##  [343,] 0.9502068
    ##  [344,] 1.0712926
    ##  [345,] 0.9162811
    ##  [346,] 1.3254131
    ##  [347,] 0.8434620
    ##  [348,] 1.0427231
    ##  [349,] 1.1858428
    ##  [350,] 0.8481673
    ##  [351,] 1.0136747
    ##  [352,] 1.1605822
    ##  [353,] 0.8116563
    ##  [354,] 0.9254515
    ##  [355,] 1.0231470
    ##  [356,] 1.2058315
    ##  [357,] 1.0278270
    ##  [358,] 1.2794941
    ##  [359,] 1.0072684
    ##  [360,] 1.0563445
    ##  [361,] 0.9611395
    ##  [362,] 1.0057800
    ##  [363,] 1.0992890
    ##  [364,] 1.0688088
    ##  [365,] 1.1968326
    ##  [366,] 0.9152601
    ##  [367,] 0.9838248
    ##  [368,] 0.8716550
    ##  [369,] 1.0398081
    ##  [370,] 0.8835869
    ##  [371,] 1.0219790
    ##  [372,] 1.1051309
    ##  [373,] 0.9588232
    ##  [374,] 1.1430581
    ##  [375,] 1.1221028
    ##  [376,] 0.8576512
    ##  [377,] 0.8373645
    ##  [378,] 1.1106056
    ##  [379,] 1.0240650
    ##  [380,] 0.8473370
    ##  [381,] 1.1093801
    ##  [382,] 1.0381448
    ##  [383,] 0.8802655
    ##  [384,] 1.0134025
    ##  [385,] 1.0241061
    ##  [386,] 0.8840747
    ##  [387,] 0.8209148
    ##  [388,] 0.9976651
    ##  [389,] 0.8821299
    ##  [390,] 1.1089728
    ##  [391,] 0.9896807
    ##  [392,] 0.8835680
    ##  [393,] 1.4451137
    ##  [394,] 0.9152601
    ##  [395,] 0.9859505
    ##  [396,] 0.7650854
    ##  [397,] 1.0935289
    ##  [398,] 0.9924519
    ##  [399,] 1.2382308
    ##  [400,] 0.9545349
    ##  [401,] 0.8999280
    ##  [402,] 0.8731724
    ##  [403,] 1.0798853
    ##  [404,] 0.8185499
    ##  [405,] 0.9770904
    ##  [406,] 0.9827484
    ##  [407,] 1.4520800
    ##  [408,] 0.9085890
    ##  [409,] 0.9916907
    ##  [410,] 0.9645999
    ##  [411,] 0.8829768
    ##  [412,] 0.7383443
    ##  [413,] 1.0016560
    ##  [414,] 0.8138443
    ##  [415,] 0.8367107
    ##  [416,] 1.2183606
    ##  [417,] 1.0139467
    ##  [418,] 0.9199714
    ##  [419,] 1.0967094
    ##  [420,] 0.7951333
    ##  [421,] 0.9540687
    ##  [422,] 1.0171745
    ##  [423,] 0.9766869
    ##  [424,] 1.2033464
    ##  [425,] 1.0111583
    ##  [426,] 1.1036381
    ##  [427,] 1.1504001
    ##  [428,] 0.8153767
    ##  [429,] 1.1373134
    ##  [430,] 0.9975559
    ##  [431,] 0.9803131
    ##  [432,] 1.0518381
    ##  [433,] 0.7540891
    ##  [434,] 1.0798853
    ##  [435,] 0.8547249
    ##  [436,] 1.2900730
    ##  [437,] 0.8284750
    ##  [438,] 1.0133659
    ##  [439,] 1.1387160
    ##  [440,] 0.9731777
    ##  [441,] 1.0028845
    ##  [442,] 0.9351232
    ##  [443,] 1.1397567
    ##  [444,] 1.0393558
    ##  [445,] 0.8585953
    ##  [446,] 1.1089728
    ##  [447,] 1.0687103
    ##  [448,] 1.0955846
    ##  [449,] 0.9788132
    ##  [450,] 1.3529435
    ##  [451,] 1.1046616
    ##  [452,] 1.2229756
    ##  [453,] 1.0289683
    ##  [454,] 1.0178981
    ##  [455,] 0.8423008
    ##  [456,] 1.0078104
    ##  [457,] 1.1416692
    ##  [458,] 1.0768972
    ##  [459,] 1.1971961
    ##  [460,] 1.2808551
    ##  [461,] 0.9551722
    ##  [462,] 1.1391317
    ##  [463,] 0.8813392
    ##  [464,] 0.9388756
    ##  [465,] 1.0209735
    ##  [466,] 0.8500055
    ##  [467,] 0.8761623
    ##  [468,] 0.8328974
    ##  [469,] 0.9910421
    ##  [470,] 0.9613304
    ##  [471,] 1.0219430
    ##  [472,] 0.8185384
    ##  [473,] 0.9922315
    ##  [474,] 0.7874424
    ##  [475,] 1.0099048
    ##  [476,] 1.0365211
    ##  [477,] 0.8012177
    ##  [478,] 0.8547249
    ##  [479,] 1.0138764
    ##  [480,] 1.0554674
    ##  [481,] 0.9465700
    ##  [482,] 0.8146654
    ##  [483,] 0.8859677
    ##  [484,] 1.1619438
    ##  [485,] 1.2176528
    ##  [486,] 0.8937220
    ##  [487,] 0.9057638
    ##  [488,] 1.1198163
    ##  [489,] 1.0649420
    ##  [490,] 0.9842195
    ##  [491,] 0.9903903
    ##  [492,] 1.0268687
    ##  [493,] 1.2675345
    ##  [494,] 0.9819008
    ##  [495,] 0.9714803
    ##  [496,] 0.8055559
    ##  [497,] 1.0449307
    ##  [498,] 1.5886831
    ##  [499,] 1.1290384
    ##  [500,] 0.9639098
    ##  [501,] 1.1223285
    ##  [502,] 1.0289683
    ##  [503,] 1.0055137
    ##  [504,] 0.8356864
    ##  [505,] 1.1605822
    ##  [506,] 0.9896807
    ##  [507,] 0.9094614
    ##  [508,] 1.0282440
    ##  [509,] 0.9385461
    ##  [510,] 1.1027348
    ##  [511,] 0.9731777
    ##  [512,] 0.8352599
    ##  [513,] 0.7060252
    ##  [514,] 0.9150637
    ##  [515,] 0.7999654
    ##  [516,] 0.9110463
    ##  [517,] 0.9663857
    ##  [518,] 1.3598460
    ##  [519,] 0.8929726
    ##  [520,] 1.2387308
    ##  [521,] 0.8120136
    ##  [522,] 1.0667912
    ##  [523,] 1.0083285
    ##  [524,] 1.3331708
    ##  [525,] 1.0491063
    ##  [526,] 0.8034662
    ##  [527,] 0.9898347
    ##  [528,] 0.9915691
    ##  [529,] 0.7705444
    ##  [530,] 1.0602969
    ##  [531,] 1.0099048
    ##  [532,] 0.8534052
    ##  [533,] 0.9366032
    ##  [534,] 0.8576512
    ##  [535,] 0.8326483
    ##  [536,] 1.0226071
    ##  [537,] 0.9530350
    ##  [538,] 1.3373899
    ##  [539,] 1.0383295
    ##  [540,] 1.1919689
    ##  [541,] 1.0405876
    ##  [542,] 1.0012430
    ##  [543,] 0.9071649
    ##  [544,] 0.9910483
    ##  [545,] 0.8425943
    ##  [546,] 1.0383740
    ##  [547,] 0.8932073
    ##  [548,] 0.8428152
    ##  [549,] 1.3617159
    ##  [550,] 0.8660866
    ##  [551,] 0.8369443
    ##  [552,] 0.8357488
    ##  [553,] 1.0377126
    ##  [554,] 0.9910421
    ##  [555,] 1.0222623
    ##  [556,] 1.0565325
    ##  [557,] 1.0673669
    ##  [558,] 0.9059748
    ##  [559,] 0.9705312
    ##  [560,] 0.9967659
    ##  [561,] 1.0241437
    ##  [562,] 1.1158413
    ##  [563,] 0.7864067
    ##  [564,] 0.9585272
    ##  [565,] 1.0251436
    ##  [566,] 1.0405907
    ##  [567,] 0.8763776
    ##  [568,] 0.8878645
    ##  [569,] 1.1593809
    ##  [570,] 1.4619103
    ##  [571,] 1.0804163
    ##  [572,] 1.0259102
    ##  [573,] 0.7831553
    ##  [574,] 0.9491554
    ##  [575,] 0.8412403
    ##  [576,] 1.3134518
    ##  [577,] 1.1010216
    ##  [578,] 1.0798853
    ##  [579,] 1.2082471
    ##  [580,] 1.1363999
    ##  [581,] 0.9248986
    ##  [582,] 0.9753959
    ##  [583,] 0.8365371
    ##  [584,] 0.8035977
    ##  [585,] 0.9236359
    ##  [586,] 0.9776911
    ##  [587,] 0.9236738
    ##  [588,] 0.9720047
    ##  [589,] 1.0663596
    ##  [590,] 0.8119721
    ##  [591,] 0.8908214
    ##  [592,] 0.9071649
    ##  [593,] 1.0240650
    ##  [594,] 1.0696500
    ##  [595,] 0.8778336
    ##  [596,] 0.9713975
    ##  [597,] 0.8712146
    ##  [598,] 1.0442704
    ##  [599,] 1.0401034
    ##  [600,] 0.8473370
    ##  [601,] 0.9313723
    ##  [602,] 0.8366617
    ##  [603,] 0.8595375
    ##  [604,] 0.8747163
    ##  [605,] 1.1769093
    ##  [606,] 1.0882170
    ##  [607,] 1.1319795
    ##  [608,] 1.0634790
    ##  [609,] 0.8585953
    ##  [610,] 1.0181088
    ##  [611,] 1.3693327
    ##  [612,] 0.9252879
    ##  [613,] 1.1858428
    ##  [614,] 0.9764786
    ##  [615,] 0.9250911
    ##  [616,] 0.8952416
    ##  [617,] 0.8826861
    ##  [618,] 0.7651444
    ##  [619,] 0.9434849
    ##  [620,] 1.0364501
    ##  [621,] 0.8874711
    ##  [622,] 0.7357947
    ##  [623,] 0.8078958
    ##  [624,] 1.0104474
    ##  [625,] 1.0506498
    ##  [626,] 0.9746421
    ##  [627,] 0.9932940
    ##  [628,] 1.0494605
    ##  [629,] 1.1223285
    ##  [630,] 1.1605822
    ##  [631,] 0.8959380
    ##  [632,] 0.9276242
    ##  [633,] 0.8576512
    ##  [634,] 1.2682253
    ##  [635,] 1.0884007
    ##  [636,] 0.7981683
    ##  [637,] 0.9191264
    ##  [638,] 1.0402614
    ##  [639,] 0.8841548
    ##  [640,] 0.7987765
    ##  [641,] 0.9071649
    ##  [642,] 1.1992392
    ##  [643,] 0.8720830
    ##  [644,] 1.0114224
    ##  [645,] 0.8797410
    ##  [646,] 0.8224442
    ##  [647,] 0.9376966
    ##  [648,] 0.9477343
    ##  [649,] 1.1281382
    ##  [650,] 0.8230913
    ##  [651,] 1.0081480
    ##  [652,] 0.9776911
    ##  [653,] 0.8926838
    ##  [654,] 1.1020888
    ##  [655,] 0.9599800
    ##  [656,] 0.8412203
    ##  [657,] 0.9282987
    ##  [658,] 1.0935289
    ##  [659,] 0.8096180
    ##  [660,] 1.1319758
    ##  [661,] 0.8171313
    ##  [662,] 0.9588232
    ##  [663,] 0.9208192
    ##  [664,] 1.0416722
    ##  [665,] 0.9388427
    ##  [666,] 1.1032447
    ##  [667,] 1.2473257
    ##  [668,] 1.0297907
    ##  [669,] 0.8369443
    ##  [670,] 0.9462332
    ##  [671,] 1.0456031
    ##  [672,] 1.0687241
    ##  [673,] 0.8967967
    ##  [674,] 1.0398462
    ##  [675,] 1.1254892
    ##  [676,] 0.9810538
    ##  [677,] 0.8185384
    ##  [678,] 0.9313723
    ##  [679,] 1.1073400
    ##  [680,] 0.9757056
    ##  [681,] 0.8222756
    ##  [682,] 0.9595837
    ##  [683,] 1.2078699
    ##  [684,] 1.0218927
    ##  [685,] 1.0402552
    ##  [686,] 0.8410361
    ##  [687,] 1.0453418
    ##  [688,] 1.0115973
    ##  [689,] 0.8907790
    ##  [690,] 1.0823235
    ##  [691,] 0.8354013
    ##  [692,] 1.1316895
    ##  [693,] 0.8343636
    ##  [694,] 0.8926433
    ##  [695,] 0.9608069
    ##  [696,] 0.9766196
    ##  [697,] 1.0611415
    ##  [698,] 1.0339571
    ##  [699,] 0.8385218
    ##  [700,] 0.8594651
    ##  [701,] 1.3322498
    ##  [702,] 1.0228481
    ##  [703,] 1.1302534
    ##  [704,] 1.0844454
    ##  [705,] 1.5506018
    ##  [706,] 1.2642162
    ##  [707,] 1.1561520
    ##  [708,] 1.1428520
    ##  [709,] 1.1440082
    ##  [710,] 1.0055137
    ##  [711,] 0.8728360
    ##  [712,] 1.0228481
    ##  [713,] 1.1197018
    ##  [714,] 0.8138443
    ##  [715,] 0.7224268
    ##  [716,] 0.9953988
    ##  [717,] 1.0053742
    ##  [718,] 1.2220610
    ##  [719,] 0.9608069
    ##  [720,] 0.9610431
    ##  [721,] 0.9592685
    ##  [722,] 1.0487984
    ##  [723,] 0.8593918
    ##  [724,] 0.8185384
    ##  [725,] 1.1561520
    ##  [726,] 0.9145883
    ##  [727,] 1.0054597
    ##  [728,] 1.1247513
    ##  [729,] 0.9764786
    ##  [730,] 0.8386403
    ##  [731,] 1.1858428
    ##  [732,] 0.9999463
    ##  [733,] 0.8588020
    ##  [734,] 0.9059748
    ##  [735,] 1.0858316
    ##  [736,] 0.9664086
    ##  [737,] 1.1412329
    ##  [738,] 0.9954901
    ##  [739,] 0.9971168
    ##  [740,] 0.7102493
    ##  [741,] 0.9743698
    ##  [742,] 0.8666307
    ##  [743,] 0.9765936
    ##  [744,] 1.0859253
    ##  [745,] 1.2038557
    ##  [746,] 0.9913084
    ##  [747,] 0.7583733
    ##  [748,] 0.9774370
    ##  [749,] 1.0684848
    ##  [750,] 0.8306064
    ##  [751,] 1.0797307
    ##  [752,] 0.9749662
    ##  [753,] 0.8547249
    ##  [754,] 1.0544406
    ##  [755,] 0.9401136
    ##  [756,] 1.3617159
    ##  [757,] 0.9169689
    ##  [758,] 1.0589150
    ##  [759,] 1.0398341
    ##  [760,] 0.7657331
    ##  [761,] 0.9095234
    ##  [762,] 0.8494926
    ##  [763,] 0.7657611
    ##  [764,] 0.9249519
    ##  [765,] 0.8829768
    ##  [766,] 1.0688406
    ##  [767,] 0.9954901
    ##  [768,] 0.9047226
    ##  [769,] 1.1123469
    ##  [770,] 1.0493231
    ##  [771,] 0.8547249
    ##  [772,] 0.8428152
    ##  [773,] 0.7882206
    ##  [774,] 1.2229756
    ##  [775,] 0.9932940
    ##  [776,] 0.9434849
    ##  [777,] 1.0123041
    ##  [778,] 1.2906741
    ##  [779,] 0.7394144
    ##  [780,] 1.1222277
    ##  [781,] 0.9528466
    ##  [782,] 1.1814642
    ##  [783,] 0.8670756
    ##  [784,] 0.9409845
    ##  [785,] 1.0288482
    ##  [786,] 1.0981311
    ##  [787,] 0.7831553
    ##  [788,] 0.8222756
    ##  [789,] 1.0944868
    ##  [790,] 1.1431534
    ##  [791,] 1.1633664
    ##  [792,] 0.8670756
    ##  [793,] 0.9749662
    ##  [794,] 1.0992890
    ##  [795,] 0.8131730
    ##  [796,] 0.8821051
    ##  [797,] 0.9589330
    ##  [798,] 0.8003510
    ##  [799,] 1.1305495
    ##  [800,] 0.8142488
    ##  [801,] 1.2519054
    ##  [802,] 0.8727798
    ##  [803,] 1.1223893
    ##  [804,] 1.0769137
    ##  [805,] 0.9161750
    ##  [806,] 1.0595208
    ##  [807,] 1.0673669
    ##  [808,] 0.7874424
    ##  [809,] 0.8802655
    ##  [810,] 1.0231470
    ##  [811,] 1.0052905
    ##  [812,] 1.1034829
    ##  [813,] 0.9149193
    ##  [814,] 1.0859253
    ##  [815,] 1.2448421
    ##  [816,] 1.0507971
    ##  [817,] 1.0530661
    ##  [818,] 0.9303378
    ##  [819,] 0.9121721
    ##  [820,] 1.0742182
    ##  [821,] 0.9288961
    ##  [822,] 1.1926923
    ##  [823,] 0.9439292
    ##  [824,] 1.0352638
    ##  [825,] 0.9513751
    ##  [826,] 1.1609070
    ##  [827,] 0.8939868
    ##  [828,] 1.0237683
    ##  [829,] 0.9944821
    ##  [830,] 1.0798853
    ##  [831,] 0.9485368
    ##  [832,] 0.9789779
    ##  [833,] 1.1093801
    ##  [834,] 0.8859430
    ##  [835,] 1.2006358
    ##  [836,] 0.9153649
    ##  [837,] 1.0123041
    ##  [838,] 0.8681251
    ##  [839,] 1.0169800
    ##  [840,] 0.8181501
    ##  [841,] 0.7691720
    ##  [842,] 1.0146761
    ##  [843,] 1.1459109
    ##  [844,] 0.9700898
    ##  [845,] 0.8898014
    ##  [846,] 1.2238573
    ##  [847,] 1.0355968
    ##  [848,] 1.0158206
    ##  [849,] 0.8859677
    ##  [850,] 0.9122456
    ##  [851,] 0.8882063
    ##  [852,] 0.9551722
    ##  [853,] 0.8897459
    ##  [854,] 0.9398469
    ##  [855,] 0.9852446
    ##  [856,] 0.9778046
    ##  [857,] 1.0393898
    ##  [858,] 0.9987543
    ##  [859,] 1.3580188
    ##  [860,] 1.1785508
    ##  [861,] 1.0738233
    ##  [862,] 1.1065117
    ##  [863,] 0.8285434
    ##  [864,] 0.9613304
    ##  [865,] 0.9202728
    ##  [866,] 0.7921230
    ##  [867,] 1.1370994
    ##  [868,] 1.1319795
    ##  [869,] 1.2887738
    ##  [870,] 1.0268119
    ##  [871,] 1.0228636
    ##  [872,] 0.9462332
    ##  [873,] 1.0592041
    ##  [874,] 0.9604953
    ##  [875,] 0.7994579
    ##  [876,] 1.1643890
    ##  [877,] 1.0427231
    ##  [878,] 1.0078104
    ##  [879,] 1.1102521
    ##  [880,] 0.9644497
    ##  [881,] 1.1396027
    ##  [882,] 0.8665142
    ##  [883,] 1.1854572
    ##  [884,] 1.5456398
    ##  [885,] 0.9311948
    ##  [886,] 0.9794400
    ##  [887,] 1.0378018
    ##  [888,] 1.0064597
    ##  [889,] 1.2329342
    ##  [890,] 1.0139467
    ##  [891,] 1.0055690
    ##  [892,] 0.8283786
    ##  [893,] 1.1223285
    ##  [894,] 0.9925005
    ##  [895,] 0.8142488
    ##  [896,] 1.2663838
    ##  [897,] 0.9623227
    ##  [898,] 0.8915469
    ##  [899,] 0.9156720
    ##  [900,] 0.9495685
    ##  [901,] 1.0820510
    ##  [902,] 0.9475504
    ##  [903,] 0.9644497
    ##  [904,] 1.1143950
    ##  [905,] 1.0554955
    ##  [906,] 1.1391317
    ##  [907,] 0.9281473
    ##  [908,] 0.8378398
    ##  [909,] 1.0454433
    ##  [910,] 1.1286916
    ##  [911,] 1.2832921
    ##  [912,] 1.1019733
    ##  [913,] 1.0469116
    ##  [914,] 0.9236089
    ##  [915,] 0.8947508
    ##  [916,] 0.9176394
    ##  [917,] 1.1410957
    ##  [918,] 1.0163980
    ##  [919,] 0.9250603
    ##  [920,] 1.0521761
    ##  [921,] 1.0398341
    ##  [922,] 1.2755469
    ##  [923,] 0.9910421
    ##  [924,] 0.8590697
    ##  [925,] 0.8595375
    ##  [926,] 1.1244386
    ##  [927,] 1.0052744
    ##  [928,] 1.1270574
    ##  [929,] 0.9954868
    ##  [930,] 1.0514969
    ##  [931,] 1.1430996
    ##  [932,] 0.9972908
    ##  [933,] 0.8882063
    ##  [934,] 1.1582040
    ##  [935,] 1.4206564
    ##  [936,] 0.8530544
    ##  [937,] 0.9482863
    ##  [938,] 1.3130098
    ##  [939,] 1.0377126
    ##  [940,] 1.0564406
    ##  [941,] 0.9924318
    ##  [942,] 1.0349919
    ##  [943,] 0.8859677
    ##  [944,] 0.9889188
    ##  [945,] 0.8318230
    ##  [946,] 1.1065117
    ##  [947,] 1.0152325
    ##  [948,] 0.9469237
    ##  [949,] 1.0544406
    ##  [950,] 1.0627765
    ##  [951,] 0.8632313
    ##  [952,] 1.1254260
    ##  [953,] 1.0449349
    ##  [954,] 0.9687935
    ##  [955,] 1.0482756
    ##  [956,] 1.0003008
    ##  [957,] 0.9780336
    ##  [958,] 1.0797307
    ##  [959,] 0.9858424
    ##  [960,] 0.8917513
    ##  [961,] 0.8410361
    ##  [962,] 1.1400401
    ##  [963,] 1.1050202
    ##  [964,] 1.0514969
    ##  [965,] 0.9482863
    ##  [966,] 0.7972923
    ##  [967,] 0.9477343
    ##  [968,] 0.9800373
    ##  [969,] 0.9580279
    ##  [970,] 0.8635325
    ##  [971,] 0.9833773
    ##  [972,] 0.8568638
    ##  [973,] 0.9208192
    ##  [974,] 0.8728360
    ##  [975,] 0.9705312
    ##  [976,] 0.9603097
    ##  [977,] 1.1619438
    ##  [978,] 1.0406153
    ##  [979,] 1.0667912
    ##  [980,] 1.1101549
    ##  [981,] 1.0806575
    ##  [982,] 0.8932073
    ##  [983,] 1.1106056
    ##  [984,] 1.1067805
    ##  [985,] 0.7972923
    ##  [986,] 0.9608069
    ##  [987,] 0.9161250
    ##  [988,] 1.1527582
    ##  [989,] 0.9953988
    ##  [990,] 1.1942718
    ##  [991,] 0.9592120
    ##  [992,] 0.9910421
    ##  [993,] 0.9258579
    ##  [994,] 1.0240650
    ##  [995,] 1.1926923
    ##  [996,] 0.9317093
    ##  [997,] 1.2256606
    ##  [998,] 0.9639098
    ##  [999,] 1.1399516
    ## 
    ## $model.matrix
    ##    (Intercept) site1 site2
    ## 1            1     1     0
    ## 2            1     1     0
    ## 3            1     1     0
    ## 4            1     1     0
    ## 5            1     0     1
    ## 6            1     0     1
    ## 7            1     0     1
    ## 8            1    -1    -1
    ## 9            1    -1    -1
    ## 10           1    -1    -1
    ## 
    ## $terms
    ## dist_tab_assay ~ site
    ## attr(,"variables")
    ## list(dist_tab_assay, site)
    ## attr(,"factors")
    ##                site
    ## dist_tab_assay    0
    ## site              1
    ## attr(,"term.labels")
    ## [1] "site"
    ## attr(,"order")
    ## [1] 1
    ## attr(,"intercept")
    ## [1] 1
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>
    ## 
    ## attr(,"class")
    ## [1] "adonis"

``` r
anova(betadisper(dist_tab_assay,samples$site))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq F value Pr(>F)
    ## Groups     2  488.84  244.42  1.2463 0.3443
    ## Residuals  7 1372.80  196.11

``` r
# Exporting results
write.csv(resOrderedDF_pv_gm, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_pv_VS_gm.csv',sep=''))
write.csv(resOrderedDF_sa_gm, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_sa_VS_gm.csv',sep=''))
write.csv(resOrderedDF_pv_sa, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_pv_VS_sa.csv',sep=''))

sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] EnhancedVolcano_1.12.0      ashr_2.2-54                
    ##  [3] apeglm_1.16.0               tximport_1.22.0            
    ##  [5] ggvenn_0.1.9                dplyr_1.0.9                
    ##  [7] vegan_2.6-2                 lattice_0.20-45            
    ##  [9] permute_0.9-7               gplots_3.1.3               
    ## [11] genefilter_1.76.0           RColorBrewer_1.1-3         
    ## [13] pheatmap_1.0.12             markdown_1.1               
    ## [15] ggrepel_0.9.1               ggplot2_3.3.6              
    ## [17] BiocManager_1.30.18         devtools_2.4.3             
    ## [19] usethis_2.1.6               DESeq2_1.34.0              
    ## [21] SummarizedExperiment_1.24.0 Biobase_2.54.0             
    ## [23] MatrixGenerics_1.6.0        matrixStats_0.62.0         
    ## [25] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
    ## [27] IRanges_2.28.0              S4Vectors_0.32.4           
    ## [29] BiocGenerics_0.40.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] plyr_1.8.7             splines_4.1.2          BiocParallel_1.28.3   
    ##   [4] digest_0.6.29          invgamma_1.1           htmltools_0.5.2       
    ##   [7] SQUAREM_2021.1         fansi_1.0.3            magrittr_2.0.3        
    ##  [10] memoise_2.0.1          cluster_2.1.3          tzdb_0.3.0            
    ##  [13] remotes_2.4.2          Biostrings_2.62.0      readr_2.1.2           
    ##  [16] annotate_1.72.0        extrafont_0.18         vroom_1.5.7           
    ##  [19] extrafontdb_1.0        bdsmatrix_1.3-6        prettyunits_1.1.1     
    ##  [22] colorspace_2.0-3       blob_1.2.3             xfun_0.31             
    ##  [25] callr_3.7.0            crayon_1.5.1           RCurl_1.98-1.6        
    ##  [28] survival_3.3-1         glue_1.6.2             gtable_0.3.0          
    ##  [31] zlibbioc_1.40.0        XVector_0.34.0         DelayedArray_0.20.0   
    ##  [34] proj4_1.0-11           pkgbuild_1.3.1         Rttf2pt1_1.3.10       
    ##  [37] maps_3.4.0             scales_1.2.0           mvtnorm_1.1-3         
    ##  [40] DBI_1.1.2              Rcpp_1.0.8.3           xtable_1.8-4          
    ##  [43] emdbook_1.3.12         bit_4.0.4              truncnorm_1.0-8       
    ##  [46] httr_1.4.3             ellipsis_0.3.2         farver_2.1.0          
    ##  [49] pkgconfig_2.0.3        XML_3.99-0.9           locfit_1.5-9.5        
    ##  [52] utf8_1.2.2             labeling_0.4.2         tidyselect_1.1.2      
    ##  [55] rlang_1.0.2            AnnotationDbi_1.56.2   munsell_0.5.0         
    ##  [58] tools_4.1.2            cachem_1.0.6           cli_3.3.0             
    ##  [61] generics_0.1.2         RSQLite_2.2.14         evaluate_0.15         
    ##  [64] stringr_1.4.0          fastmap_1.1.0          yaml_2.3.5            
    ##  [67] processx_3.5.3         knitr_1.39             bit64_4.0.5           
    ##  [70] fs_1.5.2               caTools_1.18.2         purrr_0.3.4           
    ##  [73] KEGGREST_1.34.0        nlme_3.1-157           ash_1.0-15            
    ##  [76] ggrastr_1.0.1          brio_1.1.3             compiler_4.1.2        
    ##  [79] rstudioapi_0.13        beeswarm_0.4.0         curl_4.3.2            
    ##  [82] png_0.1-7              testthat_3.1.4         tibble_3.1.7          
    ##  [85] geneplotter_1.72.0     stringi_1.7.6          ps_1.7.0              
    ##  [88] desc_1.4.1             ggalt_0.4.0            Matrix_1.4-1          
    ##  [91] vctrs_0.4.1            pillar_1.7.0           lifecycle_1.0.1       
    ##  [94] bitops_1.0-7           irlba_2.3.5            R6_2.5.1              
    ##  [97] KernSmooth_2.23-20     vipor_0.4.5            sessioninfo_1.2.2     
    ## [100] MASS_7.3-57            gtools_3.9.2.1         assertthat_0.2.1      
    ## [103] pkgload_1.2.4          rprojroot_2.0.3        withr_2.5.0           
    ## [106] GenomeInfoDbData_1.2.7 mgcv_1.8-40            parallel_4.1.2        
    ## [109] hms_1.1.1              coda_0.19-4            rmarkdown_2.14        
    ## [112] mixsqp_0.3-43          bbmle_1.0.25           numDeriv_2016.8-1.1   
    ## [115] ggbeeswarm_0.6.0
