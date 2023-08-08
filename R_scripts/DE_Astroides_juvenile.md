---
title: "DE_Astroides_adult_juveniles"
author: "Marc Meynadier"
date: "6/3/2022"
output: 
  html_document: 
    keep_md: true 
---


```r
# Differential expression on Kallisto data 

# Juveniles

# Packages and dependence
packageCheckClassic <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','pheatmap','RColorBrewer','genefilter','gplots','vegan','dplyr'))
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
#BiocManager::install('limma')
#devtools::install_github('cran/GMD')
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
```

```
## Skipping install of 'ggvenn' from a github remote, the SHA1 (5e1007f0) has not changed since last install.
##   Use `force = TRUE` to force installation
```

```r
library('ggvenn')
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
library('limma')
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
```

```
## ℹ SHA-1 hash of file is 015fc0457e61e3e93a903e69a24d96d2dac7b9fb
```

```r
# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
candidateGenes<-read.csv('candidateGenes.csv',header=T,sep=',')
samples<-read.table('tximport_design_juvenile.txt',header=T)
samples2<-read.table('tximport_design_juvenile2.txt',header=T)
samples3<-read.table('tximport_design_juvenile3.txt',header=T)
samplesNatSim<-read.table('tximport_design_juvenile_naturalSimulation.txt',header=T)
dataPath<-'/Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/juvenile/juv'
outputPath<-'/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/output/DESeq2/annotatedGenome/juvenile/juv/'
setwd(dataPath)
data<-list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(data, read.table, skip = 4)
raw_counts <- as.data.frame(sapply(counts.files, function(x) x[ , 4]))
data <- gsub( "Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/juvenile/juv", "", data )
data <- gsub( "_ReadsPerGene.out.tab", "", data )
data <- gsub( "./", "", data )
colnames(raw_counts) <- data
row.names(raw_counts) <- counts.files[[1]]$V1

raw_counts_2 <- raw_counts[,grep("gm_amb|gm_low|sp_amb|sp_low", colnames(raw_counts))] 
raw_counts_natSim <- raw_counts[,grep("gm_low|sp_amb", colnames(raw_counts))] 

# DDS object 

dds<-DESeqDataSetFromMatrix(countData = raw_counts,colData=samples,design= ~site + pH)
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors
```

```r
dds2<-DESeqDataSetFromMatrix(countData = raw_counts_2,colData=samples2,design= ~site + pH)
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors
```

```r
dds3<-DESeqDataSetFromMatrix(countData = raw_counts,colData=samples3,design= ~site_pH)
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors
```

```r
ddsNatSim<-DESeqDataSetFromMatrix(countData = raw_counts_natSim,colData=samplesNatSim,design= ~site)
```

```
## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors
```

```r
# If data importation using Kallisto
#tx2gene<-read.table('tx2gene_adultTranscriptome',header=T)
#scriptPath <- sub("/[^/]+$", "", scriptPath)
#scriptPath <- sub("/[^/]+$", "", scriptPath)
#dataPath<-'/data/net/6_kallisto/adultTranscriptome/juvenile'
#outputPath<-paste(scriptPath,'/output/DESeq2/adultTranscriptome/juvenile/',sep='')
#wdPath<-paste(scriptPath,dataPath,sep='')
#setwd(wdPath)
#files<-paste0(samples$samples,'.tsv')
#files2<-paste0(samples2$samples,'.tsv')
#files3<-paste0(samples3$samples,'.tsv')
#filesNatSim<-paste0(samplesNatSim$samples,'.tsv')
#names(files)<-samples$samples
#names(files2)<-samples2$samples
#names(files3)<-samples3$samples
#names(filesNatSim)<-samplesNatSim$samples
#txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
#txi2<-tximport(files = files2,type='kallisto',tx2gene = tx2gene)
#txi3<-tximport(files = files3,type='kallisto',tx2gene = tx2gene)
#txiNatSim<-tximport(files = filesNatSim,type='kallisto',tx2gene = tx2gene)
#names(txi)
#names(txi2)
#names(txi3)
#names(txiNatSim)
#head(txi$counts)
#head(txi2$counts)
#head(txi3$counts)
#head(txiNatSim$counts)
#dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site + pH)
#dds2<-DESeqDataSetFromTximport(txi2,colData=samples2,design= ~site + pH)
#dds3<-DESeqDataSetFromTximport(txi3,colData=samples3,design= ~site_pH)
#ddsNatSim<-DESeqDataSetFromTximport(txiNatSim,colData=samplesNatSim,design= ~Site)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]
keep <- rowSums(counts(dds2)) >= 10 
dds2 <- dds2[keep,]
keep <- rowSums(counts(dds3)) >= 10 
dds3 <- dds3[keep,]
keep <- rowSums(counts(ddsNatSim)) >= 10 
ddsNatSim <- ddsNatSim[keep,]

# Differential expression analysis
dds<-DESeq(dds)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```r
dds2<-DESeq(dds2)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```r
dds3<-DESeq(dds3)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```r
ddsNatSim<-DESeq(ddsNatSim)
```

```
## estimating size factors
```

```
## estimating dispersions
```

```
## gene-wise dispersion estimates
```

```
## mean-dispersion relationship
```

```
## final dispersion estimates
```

```
## fitting model and testing
```

```r
cbind(resultsNames(dds))
```

```
##      [,1]                       
## [1,] "Intercept"                
## [2,] "site_SP_vs_GM"            
## [3,] "pH_extreme_low_vs_ambient"
## [4,] "pH_low_vs_ambient"
```

```r
cbind(resultsNames(dds2))
```

```
##      [,1]               
## [1,] "Intercept"        
## [2,] "site_SP_vs_GM"    
## [3,] "pH_low_vs_ambient"
```

```r
cbind(resultsNames(dds3))
```

```
##      [,1]                                  
## [1,] "Intercept"                           
## [2,] "site_pH_GM_extreme_low_vs_GM_ambient"
## [3,] "site_pH_GM_low_vs_GM_ambient"        
## [4,] "site_pH_SP_ambient_vs_GM_ambient"    
## [5,] "site_pH_SP_low_vs_GM_ambient"
```

```r
cbind(resultsNames(ddsNatSim))
```

```
##      [,1]           
## [1,] "Intercept"    
## [2,] "site_SP_vs_GM"
```

```r
sp_VS_gm_global<-results(dds, contrast=c("site","SP","GM"), alpha = 0.05)
sp_VS_gm<-results(dds2, contrast=c("site","SP","GM"), alpha = 0.05)
amb_VS_ext<-results(dds, contrast=c("pH","ambient","extreme_low"), alpha = 0.05)
amb_VS_low<-results(dds, contrast=c("pH","ambient","low"), alpha = 0.05)
low_VS_ext<-results(dds, contrast=c("pH","low","extreme_low"), alpha = 0.05)
gm_amb_VS_gm_low<-results(dds3, contrast=c("site_pH","GM_ambient","GM_low"), alpha = 0.05)
gm_amb_VS_gm_ext<-results(dds3, contrast=c("site_pH","GM_ambient","GM_extreme_low"), alpha = 0.05)
gm_low_VS_gm_ext<-results(dds3, contrast=c("site_pH","GM_ambient","GM_extreme_low"), alpha = 0.05)
sp_amb_VS_sp_low<-results(dds3, contrast=c("site_pH","SP_ambient","SP_low"), alpha = 0.05)
gm_amb_VS_sp_amb<-results(dds3, contrast=c("site_pH","GM_ambient","SP_ambient"), alpha = 0.05)
gm_low_VS_sp_low<-results(dds3, contrast=c("site_pH","GM_low","SP_low"), alpha = 0.05)
sp_amb_VS_gm_low_natSim<-results(ddsNatSim,contrast=c("site","GM","SP"),alpha = 0.05)
summary(sp_VS_gm_global)
```

```
## 
## out of 14704 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 40, 0.27%
## LFC < 0 (down)     : 42, 0.29%
## outliers [1]       : 205, 1.4%
## low counts [2]     : 7084, 48%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
summary(sp_VS_gm)
```

```
## 
## out of 13516 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 49, 0.36%
## LFC < 0 (down)     : 43, 0.32%
## outliers [1]       : 226, 1.7%
## low counts [2]     : 6751, 50%
## (mean count < 2)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
summary(amb_VS_ext)
```

```
## 
## out of 14704 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 3, 0.02%
## LFC < 0 (down)     : 0, 0%
## outliers [1]       : 205, 1.4%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
summary(amb_VS_low)
```

```
## 
## out of 14704 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 3, 0.02%
## LFC < 0 (down)     : 1, 0.0068%
## outliers [1]       : 205, 1.4%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
summary(low_VS_ext)
```

```
## 
## out of 14704 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 3, 0.02%
## LFC < 0 (down)     : 2, 0.014%
## outliers [1]       : 205, 1.4%
## low counts [2]     : 0, 0%
## (mean count < 0)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
summary(sp_amb_VS_gm_low_natSim)
```

```
## 
## out of 9771 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 14, 0.14%
## LFC < 0 (down)     : 40, 0.41%
## outliers [1]       : 137, 1.4%
## low counts [2]     : 4545, 47%
## (mean count < 3)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
```

```r
# Exploring the results

# Results sp VS gm

#MA-plot
DESeq2::plotMA(sp_VS_gm,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nsp_VS_gm")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

```r
# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
EnhancedVolcano(data.frame(sp_VS_gm), lab = rownames(data.frame(sp_VS_gm)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_VS_gm), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-2.png)

```r
# Results ext VS amb

#MA-plot
DESeq2::plotMA(amb_VS_ext,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\namb_VS_ext")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-3.png)

```r
# Volcano plot
EnhancedVolcano(data.frame(amb_VS_ext), lab = rownames(data.frame(amb_VS_ext)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between ext and amb",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(amb_VS_ext), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-4.png)

```r
# Results low VS amb

#MA-plot
DESeq2::plotMA(amb_VS_low,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\namb_VS_low")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-5.png)

```r
# Volcano plot
EnhancedVolcano(data.frame(amb_VS_low), lab = rownames(data.frame(amb_VS_low)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between low and amb",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(amb_VS_low), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-6.png)

```r
# Results low VS ext

#MA-plot
DESeq2::plotMA(low_VS_ext,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nlow_VS_ext")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-7.png)

```r
# Volcano plot
EnhancedVolcano(data.frame(low_VS_ext), lab = rownames(data.frame(low_VS_ext)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between low and ext",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(low_VS_ext), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-8.png)

```r
# Results natural simulation

#MA-plot
DESeq2::plotMA(sp_amb_VS_gm_low_natSim,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nsp_amb_VS_gm_low")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-9.png)

```r
# Volcano plot
EnhancedVolcano(data.frame(sp_amb_VS_gm_low_natSim), lab = rownames(data.frame(sp_amb_VS_gm_low_natSim)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp_amb and gm_low",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_amb_VS_gm_low_natSim), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-10.png)

```r
# Principal Component Analysis

vsd3 = vst(dds3,blind=T)

pcaData = plotPCA(vsd3, intgroup="site_pH", 
                  returnData=TRUE,ntop=1000)
percentVar = round(100 * attr(pcaData, "percentVar"))

pcaData$site_pH = factor(pcaData$site_pH, levels=c("GM_extreme_low","GM_low","GM_ambient","SP_low","SP_ambient"))

ggplot(pcaData, aes(PC1, PC2, fill = site_pH)) + 
  geom_point(color="black",pch=21, size=5) + theme_bw() +
  scale_fill_manual(values = c("#D55E00","#E69F00","#FFFF00","#0072B2","#56B4E9")) +
  #ggtitle("Principal Component Analysis of adult corals", subtitle = "may2018 dataset") +
  theme(text = element_text(size=14), legend.position = 'bottom') +
  theme(legend.title=element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-11.png)

```r
vsdNatSim = vst(ddsNatSim,blind=T)

pcaData = plotPCA(vsdNatSim, intgroup="site", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of juvenile corals", subtitle = "Juvenile dataset - Natural conditions simulation") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-12.png)

```r
# Venn diagramm 
resOrdered_sp_VS_gm <- sp_VS_gm_global[order(sp_VS_gm_global$padj),]
resOrderedDF_sp_VS_gm <- as.data.frame(resOrdered_sp_VS_gm)
resOrderedDF_sp_VS_gm_venn <- filter(resOrderedDF_sp_VS_gm,padj < 0.05)
resOrderedDF_sp_VS_gm_venn <- list(rownames(resOrderedDF_sp_VS_gm_venn))
resOrderedDF_sp_VS_gm_venn <- unlist(resOrderedDF_sp_VS_gm_venn)

resOrdered_amb_VS_ext <- amb_VS_ext[order(amb_VS_ext$padj),]
resOrderedDF_amb_VS_ext <- as.data.frame(resOrdered_amb_VS_ext)
resOrderedDF_amb_VS_ext_venn <- filter(resOrderedDF_amb_VS_ext,padj < 0.05)
resOrderedDF_amb_VS_ext_venn <- list(rownames(resOrderedDF_amb_VS_ext_venn))
resOrderedDF_amb_VS_ext_venn <- unlist(resOrderedDF_amb_VS_ext_venn)

resOrdered_amb_VS_low <- amb_VS_low[order(amb_VS_low$padj),]
resOrderedDF_amb_VS_low <- as.data.frame(resOrdered_amb_VS_low)
resOrderedDF_amb_VS_low_venn <- filter(resOrderedDF_amb_VS_low,padj < 0.05)
resOrderedDF_amb_VS_low_venn <- list(rownames(resOrderedDF_amb_VS_low_venn))
resOrderedDF_amb_VS_low_venn <- unlist(resOrderedDF_amb_VS_low_venn)

resOrdered_low_VS_ext <- low_VS_ext[order(low_VS_ext$padj),]
resOrderedDF_low_VS_ext <- as.data.frame(resOrdered_low_VS_ext)
resOrderedDF_low_VS_ext_venn <- filter(resOrderedDF_low_VS_ext,padj < 0.05)
resOrderedDF_low_VS_ext_venn <- list(rownames(resOrderedDF_low_VS_ext_venn))
resOrderedDF_low_VS_ext_venn <- unlist(resOrderedDF_low_VS_ext_venn)

x = list('sp VS gm' = resOrderedDF_sp_VS_gm_venn, 'amb VS ext' = resOrderedDF_amb_VS_ext_venn,
         'amb VS low' = resOrderedDF_amb_VS_low_venn, 'low VS ext' = resOrderedDF_low_VS_ext_venn)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#009E73"),
  stroke_size = 0.4, set_name_size = 4
)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-13.png)

```r
resOrdered_sp_VS_gm <- sp_VS_gm[order(sp_VS_gm$padj),]
resOrderedDF_sp_VS_gm <- as.data.frame(resOrdered_sp_VS_gm)
resOrderedDF_sp_VS_gm_venn <- filter(resOrderedDF_sp_VS_gm,padj < 0.05)
resOrderedDF_sp_VS_gm_venn <- list(rownames(resOrderedDF_sp_VS_gm_venn))
resOrderedDF_sp_VS_gm_venn <- unlist(resOrderedDF_sp_VS_gm_venn)

x = list('SP VS GM' = resOrderedDF_sp_VS_gm_venn,'ambient VS low' = resOrderedDF_amb_VS_low_venn)

ggvenn(
  x, 
  fill_color = c("#8E8E8E", "#E0E0E0"),
  stroke_size = 0.4, set_name_size = 4
)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-14.png)

```r
resOrdered_gm_amb_VS_gm_low <- gm_amb_VS_gm_low[order(gm_amb_VS_gm_low$padj),]
resOrderedDF_gm_amb_VS_gm_low <- as.data.frame(resOrdered_gm_amb_VS_gm_low)
resOrderedDF_gm_amb_VS_gm_low_venn <- filter(resOrderedDF_gm_amb_VS_gm_low,padj < 0.05)
resOrderedDF_gm_amb_VS_gm_low_venn <- list(rownames(resOrderedDF_gm_amb_VS_gm_low_venn))
resOrderedDF_gm_amb_VS_gm_low_venn <- unlist(resOrderedDF_gm_amb_VS_gm_low_venn)

resOrdered_gm_amb_VS_gm_ext <- gm_amb_VS_gm_ext[order(gm_amb_VS_gm_ext$padj),]
resOrderedDF_gm_amb_VS_gm_ext <- as.data.frame(resOrdered_gm_amb_VS_gm_ext)
resOrderedDF_gm_amb_VS_gm_ext_venn <- filter(resOrderedDF_gm_amb_VS_gm_ext,padj < 0.05)
resOrderedDF_gm_amb_VS_gm_ext_venn <- list(rownames(resOrderedDF_gm_amb_VS_gm_ext_venn))
resOrderedDF_gm_amb_VS_gm_ext_venn <- unlist(resOrderedDF_gm_amb_VS_gm_ext_venn)

resOrdered_gm_low_VS_gm_ext <- gm_low_VS_gm_ext[order(gm_low_VS_gm_ext$padj),]
resOrderedDF_gm_low_VS_gm_ext <- as.data.frame(resOrdered_gm_low_VS_gm_ext)
resOrderedDF_gm_low_VS_gm_ext_venn <- filter(resOrderedDF_gm_low_VS_gm_ext,padj < 0.05)
resOrderedDF_gm_low_VS_gm_ext_venn <- list(rownames(resOrderedDF_gm_low_VS_gm_ext_venn))
resOrderedDF_gm_low_VS_gm_ext_venn <- unlist(resOrderedDF_gm_low_VS_gm_ext_venn)

resOrdered_sp_amb_VS_sp_low <- sp_amb_VS_sp_low[order(sp_amb_VS_sp_low$padj),]
resOrderedDF_sp_amb_VS_sp_low <- as.data.frame(resOrdered_sp_amb_VS_sp_low)
resOrderedDF_sp_amb_VS_sp_low_venn <- filter(resOrderedDF_sp_amb_VS_sp_low,padj < 0.05)
resOrderedDF_sp_amb_VS_sp_low_venn <- list(rownames(resOrderedDF_sp_amb_VS_sp_low_venn))
resOrderedDF_sp_amb_VS_sp_low_venn <- unlist(resOrderedDF_sp_amb_VS_sp_low_venn)

resOrdered_gm_low_VS_sp_low <- gm_low_VS_sp_low[order(gm_low_VS_sp_low$padj),]
resOrderedDF_gm_low_VS_sp_low <- as.data.frame(resOrdered_gm_low_VS_sp_low)
resOrderedDF_gm_low_VS_sp_low_venn <- filter(resOrderedDF_gm_low_VS_sp_low,padj < 0.05)
resOrderedDF_gm_low_VS_sp_low_venn <- list(rownames(resOrderedDF_gm_low_VS_sp_low_venn))
resOrderedDF_gm_low_VS_sp_low_venn <- unlist(resOrderedDF_gm_low_VS_sp_low_venn)

resOrdered_gm_amb_VS_sp_amb <- gm_amb_VS_sp_amb[order(gm_amb_VS_sp_amb$padj),]
resOrderedDF_gm_amb_VS_sp_amb <- as.data.frame(resOrdered_gm_amb_VS_sp_amb)
resOrderedDF_gm_amb_VS_sp_amb_venn <- filter(resOrderedDF_gm_amb_VS_sp_amb,padj < 0.05)
resOrderedDF_gm_amb_VS_sp_amb_venn <- list(rownames(resOrderedDF_gm_amb_VS_sp_amb_venn))
resOrderedDF_gm_amb_VS_sp_amb_venn <- unlist(resOrderedDF_gm_amb_VS_sp_amb_venn)

x = list('GM ambient VS GM low' = resOrderedDF_gm_amb_VS_gm_low_venn,'SP ambient VS SP low' = resOrderedDF_sp_amb_VS_sp_low_venn)

ggvenn(
  x, 
  fill_color = c("#EE4000", "#5CACEE"),
  stroke_size = 0.4, set_name_size = 4
)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-15.png)

```r
x = list('GM ambient VS SP ambient' = resOrderedDF_gm_amb_VS_sp_amb_venn,'GM low VS SP low' = resOrderedDF_gm_low_VS_sp_low_venn)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-16.png)

```r
x = list('GM ambient VS GM extreme low' = resOrderedDF_gm_amb_VS_gm_ext_venn,'GM low VS GM extreme low' = resOrderedDF_gm_low_VS_gm_ext_venn)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-17.png)

```r
# Candidate genes heatmap

# Global
listGenes <- candidateGenes$genes
listGenes <- gsub("^([^.]*.[^.]*).*$", "\\1", listGenes)
listGenes <- unique(listGenes)
listGenes2 <- which(rownames(vsd3) %in% listGenes)
index <- which(listGenes %in% rownames(vsd3))
candidateGenes2 <- candidateGenes[index, ] 
listProt <- candidateGenes2$pfam_annotation
listGenes3 <- candidateGenes2$genes
listGenes3 <- gsub("^([^.]*.[^.]*).*$", "\\1", listGenes3)

removeGenes = c("STRG.36428","STRG.36427","STRG.424","STRG.34997","STRG.14728","STRG.11250","STRG.37944","STRG.35059")
indexRemoveGenes = which(listGenes3 %in% removeGenes)
listProt <- listProt[-indexRemoveGenes]
listGenes3 <- listGenes3[! listGenes3 %in% c("STRG.36428","STRG.36427","STRG.424","STRG.34997","STRG.14728","STRG.11250","STRG.37944","STRG.35059")]

vsdCandidate <- vsd3[listGenes3, ]

labColName <- samples3$site_pH
colnames(vsdCandidate) <- labColName
rownames(vsdCandidate) <- paste(listProt,listGenes3,sep=" - ")

topVarGenesVsd <- order(rowVars(assay(vsdCandidate)), decreasing=TRUE)
assayVsdCandidate<-unique(assay(vsdCandidate))
pheatmap(assayVsdCandidate,scale = "row")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-18.png)

```r
# Natural conditions simulation
listGenes <- candidateGenes$genes
listGenes <- gsub("^([^.]*.[^.]*).*$", "\\1", listGenes)
listGenes <- unique(listGenes)
listGenes2 <- which(rownames(vsd) %in% listGenes)
```

```
## Error in h(simpleError(msg, call)): erreur d'ï¿½valuation de l'argument 'x' lors de la sï¿½lection d'une mï¿½thode pour la fonction 'which' : erreur d'ï¿½valuation de l'argument 'x' lors de la sï¿½lection d'une mï¿½thode pour la fonction '%in%' : erreur d'ï¿½valuation de l'argument 'x' lors de la sï¿½lection d'une mï¿½thode pour la fonction 'rownames' : objet 'vsd' introuvable
```

```r
index <- which(listGenes %in% rownames(vsdNatSim))
candidateGenes2 <- candidateGenes[index, ] 
listProt <- candidateGenes2$pfam_annotation
listGenes3 <- candidateGenes2$genes
listGenes3 <- gsub("^([^.]*.[^.]*).*$", "\\1", listGenes3)

removeGenes = c("STRG.13175","STRG.36428","STRG.424","STRG.15785","STRG.36427","STRG.34997","STRG.37944","STRG.14728","STRG.35059")
indexRemoveGenes = which(listGenes3 %in% removeGenes)
listProt <- listProt[-indexRemoveGenes]
listGenes3 <- listGenes3[! listGenes3 %in% c("STRG.13175","STRG.36428","STRG.424","STRG.15785","STRG.36427","STRG.34997","STRG.37944","STRG.14728","STRG.35059")]

vsdCandidate <- vsdNatSim[listGenes3, ]

labColName <- samples$site
colnames(vsdCandidate) <- labColName
```

```
## Error in `rownames<-`(`*tmp*`, value = value[[2]]): invalid rownames length
```

```r
rownames(vsdCandidate) <- paste(listProt,listGenes3,sep=" - ")

topVarGenesVsd <- order(rowVars(assay(vsdCandidate)), decreasing=TRUE)
assayVsdCandidate<-unique(assay(vsdCandidate))
pheatmap(assayVsdCandidate,scale = "row")
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-19.png)

```r
# Inferences statistics

count_tab_assay <- assay(vsd)
```

```
## Error in h(simpleError(msg, call)): erreur d'ï¿½valuation de l'argument 'x' lors de la sï¿½lection d'une mï¿½thode pour la fonction 'assay' : objet 'vsd' introuvable
```

```r
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site + pH, method="euclidian")
```

```
## 'adonis' will be deprecated: use 'adonis2' instead
```

```
## Error in G * t(hat): tableaux de tailles inadéquates
```

```r
anova(betadisper(dist_tab_assay,samples$site))
```

```
## Error in pts[groups == i, , drop = FALSE]: (subscript) indice logique trop long
```

```r
anova(betadisper(dist_tab_assay,samples$pH))
```

```
## Error in pts[groups == i, , drop = FALSE]: (subscript) indice logique trop long
```

```r
count_tab_assay <- assay(vsd3)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples3,dist_tab_assay ~ site_pH, method="euclidian")
```

```
## 'adonis' will be deprecated: use 'adonis2' instead
```

```
## $aov.tab
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
## site_pH    4     52252   13063  1.0757 0.22291  0.102
## Residuals 15    182157   12144         0.77709       
## Total     19    234408                 1.00000       
## 
## $call
## adonis(formula = dist_tab_assay ~ site_pH, data = samples3, method = "euclidian")
## 
## $coefficients
## NULL
## 
## $coef.sites
##                    [,1]       [,2]        [,3]        [,4]        [,5]        [,6]        [,7]      [,8]        [,9]      [,10]      [,11]      [,12]       [,13]
## (Intercept) 130.7588283 134.807073 135.8949727 134.7317003 134.2038253 133.4007881 132.2741800 280.25813 134.0538322 134.542738 139.308286 133.584563 136.2183344
## site_pH1    -32.7959082  -1.893937  -0.7900523   2.7097770   1.4986823 -32.8204387  -0.8041956  21.68713   0.3892139   1.949881 -37.526880  -4.816427  -1.9571189
## site_pH2     -5.7457333 -36.586925  -6.2378775   0.6769563  -0.6941418  -4.4511453 -35.7846603  11.87125  -1.0643103   0.586758  -2.943672 -35.624813  -3.5356458
## site_pH3     38.8363085  36.813520   4.0301841  43.2871940  39.6180995  37.5616934  36.5400653 -57.08348  42.6901600  43.688309  40.914347  38.467168   5.9447366
## site_pH4     -0.7802774   1.078353   2.3434225 -39.6673349  -4.9258506  -0.6349917  -0.4696334  12.69024 -37.3155013  -9.665478  -1.046136   1.012269  -0.2393224
##                   [,14]       [,15]       [,16]        [,17]       [,18]      [,19]       [,20]
## (Intercept) 131.0420753 132.8442956 141.2612811 134.39937668 135.6900188 132.760191 134.4122406
## site_pH1      2.9212218   2.8033797 -37.5167533   0.09090839  -4.2270862   3.887627   4.4101712
## site_pH2      0.6252415  -0.6679164  -3.9470922 -35.48490694  -4.8838737   1.774029   0.9977088
## site_pH3     43.4357915  39.5822925  40.7327886  38.38936376   7.4920157  41.695197  41.2247894
## site_pH4    -39.1315617  -6.4204502   0.2278606  -2.08646553   0.8403476 -39.488516  -9.8287073
## 
## $f.perms
##              [,1]
##    [1,] 1.0119038
##    [2,] 1.0733551
##    [3,] 0.9653052
##    [4,] 0.9595464
##    [5,] 1.0378016
##    [6,] 1.0003895
##    [7,] 1.1026315
##    [8,] 0.9712918
##    [9,] 1.0209961
##   [10,] 1.0069364
##   [11,] 1.0238364
##   [12,] 1.0379162
##   [13,] 1.0253871
##   [14,] 1.0045010
##   [15,] 0.9711923
##   [16,] 0.9708506
##   [17,] 1.0471469
##   [18,] 0.9578590
##   [19,] 1.0624692
##   [20,] 0.8938131
##   [21,] 1.0303512
##   [22,] 0.9459239
##   [23,] 1.0964458
##   [24,] 0.9811509
##   [25,] 0.9789118
##   [26,] 1.0134354
##   [27,] 0.9791706
##   [28,] 0.9788809
##   [29,] 1.0573556
##   [30,] 0.9645708
##   [31,] 0.9691987
##   [32,] 0.9352132
##   [33,] 1.0398918
##   [34,] 1.0027856
##   [35,] 1.0092266
##   [36,] 0.9474096
##   [37,] 1.0201682
##   [38,] 0.9948796
##   [39,] 1.0151470
##   [40,] 0.9718135
##   [41,] 1.0203676
##   [42,] 0.9419270
##   [43,] 1.0039612
##   [44,] 1.0066225
##   [45,] 1.0009504
##   [46,] 0.9835804
##   [47,] 1.0550251
##   [48,] 0.9718624
##   [49,] 0.9821675
##   [50,] 1.0328373
##   [51,] 1.0428273
##   [52,] 1.0818913
##   [53,] 0.9107315
##   [54,] 0.9538153
##   [55,] 1.1969454
##   [56,] 1.0875063
##   [57,] 0.9902606
##   [58,] 0.9411022
##   [59,] 0.9941504
##   [60,] 0.9501685
##   [61,] 1.0008341
##   [62,] 1.0911135
##   [63,] 0.9841305
##   [64,] 1.0710707
##   [65,] 1.1013356
##   [66,] 0.9612097
##   [67,] 0.9312107
##   [68,] 0.9300269
##   [69,] 0.9828719
##   [70,] 0.9919855
##   [71,] 1.0005168
##   [72,] 1.0082713
##   [73,] 0.9484636
##   [74,] 1.0415439
##   [75,] 1.0183914
##   [76,] 1.0661610
##   [77,] 0.9451114
##   [78,] 0.9707418
##   [79,] 1.0563331
##   [80,] 0.9555354
##   [81,] 0.9610823
##   [82,] 1.0730840
##   [83,] 1.0034462
##   [84,] 1.0024690
##   [85,] 0.9408961
##   [86,] 1.0588294
##   [87,] 0.9782302
##   [88,] 1.0433137
##   [89,] 0.9460799
##   [90,] 0.9386273
##   [91,] 0.9407847
##   [92,] 0.9482270
##   [93,] 1.0376538
##   [94,] 1.0420796
##   [95,] 1.0114705
##   [96,] 1.0945737
##   [97,] 0.9382702
##   [98,] 1.0590311
##   [99,] 0.9775514
##  [100,] 1.0797075
##  [101,] 0.9076173
##  [102,] 1.0044587
##  [103,] 1.1133191
##  [104,] 1.0134253
##  [105,] 0.9493872
##  [106,] 0.9655481
##  [107,] 0.9793494
##  [108,] 1.0760339
##  [109,] 1.1025720
##  [110,] 0.9643384
##  [111,] 0.9812822
##  [112,] 1.0281804
##  [113,] 0.9679015
##  [114,] 0.9587804
##  [115,] 1.0047270
##  [116,] 0.9834914
##  [117,] 1.0574341
##  [118,] 0.9890895
##  [119,] 1.0369425
##  [120,] 0.9190309
##  [121,] 1.0386848
##  [122,] 0.9402116
##  [123,] 1.0055099
##  [124,] 0.9248424
##  [125,] 0.9581758
##  [126,] 1.1116745
##  [127,] 0.9746253
##  [128,] 0.9592259
##  [129,] 0.9534530
##  [130,] 1.0402338
##  [131,] 1.0831652
##  [132,] 0.9983164
##  [133,] 0.9798952
##  [134,] 0.9816641
##  [135,] 1.0470933
##  [136,] 1.0444573
##  [137,] 0.9794144
##  [138,] 0.9469518
##  [139,] 1.0745768
##  [140,] 1.0257819
##  [141,] 1.0118337
##  [142,] 0.9954962
##  [143,] 0.9771947
##  [144,] 1.0416655
##  [145,] 0.9522240
##  [146,] 0.9171440
##  [147,] 0.9513210
##  [148,] 0.9275813
##  [149,] 0.9684493
##  [150,] 0.9628388
##  [151,] 0.9432265
##  [152,] 1.0159039
##  [153,] 0.9925314
##  [154,] 0.9871620
##  [155,] 0.9856287
##  [156,] 1.0080853
##  [157,] 1.0217515
##  [158,] 1.0116021
##  [159,] 1.0079220
##  [160,] 0.9771577
##  [161,] 0.9995652
##  [162,] 1.0111475
##  [163,] 0.9318303
##  [164,] 1.0745437
##  [165,] 0.9858114
##  [166,] 0.9658995
##  [167,] 1.0855430
##  [168,] 1.0387605
##  [169,] 1.1222110
##  [170,] 0.9665375
##  [171,] 1.0423052
##  [172,] 1.0744324
##  [173,] 0.9622611
##  [174,] 0.9516058
##  [175,] 1.0181889
##  [176,] 0.9403000
##  [177,] 1.0139371
##  [178,] 0.9575755
##  [179,] 1.0017725
##  [180,] 1.0268851
##  [181,] 1.0084757
##  [182,] 0.9629640
##  [183,] 1.0382747
##  [184,] 1.0443721
##  [185,] 1.0217973
##  [186,] 0.9971409
##  [187,] 0.9752433
##  [188,] 1.0546918
##  [189,] 1.0784088
##  [190,] 1.0848093
##  [191,] 1.0364977
##  [192,] 0.9519493
##  [193,] 1.0057903
##  [194,] 1.0163256
##  [195,] 0.9410364
##  [196,] 1.0616000
##  [197,] 1.0139112
##  [198,] 1.0334883
##  [199,] 1.1003879
##  [200,] 1.0391843
##  [201,] 0.9929313
##  [202,] 0.9461922
##  [203,] 1.1487608
##  [204,] 1.0194731
##  [205,] 1.0798654
##  [206,] 1.0377649
##  [207,] 0.9958725
##  [208,] 0.9558866
##  [209,] 0.9193604
##  [210,] 0.9794808
##  [211,] 1.0226000
##  [212,] 1.0036630
##  [213,] 0.9833339
##  [214,] 0.9474585
##  [215,] 0.9871509
##  [216,] 1.0438305
##  [217,] 0.9975159
##  [218,] 0.9715194
##  [219,] 0.9858694
##  [220,] 0.9844232
##  [221,] 1.0014217
##  [222,] 0.9872991
##  [223,] 1.0209478
##  [224,] 1.0597959
##  [225,] 1.1277696
##  [226,] 0.9967341
##  [227,] 1.0559309
##  [228,] 1.0012150
##  [229,] 0.8968603
##  [230,] 0.9868317
##  [231,] 1.0338769
##  [232,] 1.0170275
##  [233,] 0.9699453
##  [234,] 0.9556481
##  [235,] 1.0204175
##  [236,] 0.9289325
##  [237,] 0.9168153
##  [238,] 1.1528059
##  [239,] 0.9888308
##  [240,] 0.9772270
##  [241,] 0.9769182
##  [242,] 1.1051425
##  [243,] 1.0766953
##  [244,] 1.0407660
##  [245,] 0.9095353
##  [246,] 0.9755837
##  [247,] 0.9944533
##  [248,] 1.0104647
##  [249,] 0.9681691
##  [250,] 0.9849302
##  [251,] 0.9775541
##  [252,] 0.9950278
##  [253,] 1.0285688
##  [254,] 0.9714036
##  [255,] 1.0060756
##  [256,] 1.0346437
##  [257,] 0.9172865
##  [258,] 1.1079302
##  [259,] 1.0531788
##  [260,] 0.9999991
##  [261,] 1.0313019
##  [262,] 0.9758234
##  [263,] 1.0079053
##  [264,] 0.9604465
##  [265,] 0.9711630
##  [266,] 0.8986290
##  [267,] 1.0157010
##  [268,] 1.0605307
##  [269,] 0.9933728
##  [270,] 1.0186552
##  [271,] 1.0046114
##  [272,] 0.9921202
##  [273,] 1.0703708
##  [274,] 1.0608959
##  [275,] 1.0387850
##  [276,] 1.0169254
##  [277,] 1.1136338
##  [278,] 1.0058699
##  [279,] 1.0874683
##  [280,] 0.9481672
##  [281,] 0.9851891
##  [282,] 1.1356872
##  [283,] 1.0451481
##  [284,] 0.9715547
##  [285,] 1.0644457
##  [286,] 1.0001292
##  [287,] 0.9905541
##  [288,] 1.0303042
##  [289,] 0.9829443
##  [290,] 0.9131186
##  [291,] 0.9956182
##  [292,] 1.0328715
##  [293,] 0.9780014
##  [294,] 1.0112985
##  [295,] 0.9288326
##  [296,] 1.0167545
##  [297,] 1.0684642
##  [298,] 1.0400041
##  [299,] 0.9170360
##  [300,] 0.9698622
##  [301,] 0.9926524
##  [302,] 1.0839575
##  [303,] 0.9584320
##  [304,] 0.9324545
##  [305,] 1.0174276
##  [306,] 1.0626768
##  [307,] 1.0041821
##  [308,] 1.0813439
##  [309,] 1.0862313
##  [310,] 0.9789312
##  [311,] 0.9921255
##  [312,] 1.0087649
##  [313,] 1.1116457
##  [314,] 0.9763193
##  [315,] 1.0237799
##  [316,] 0.9717345
##  [317,] 1.0915978
##  [318,] 0.9018533
##  [319,] 0.9536722
##  [320,] 0.9798313
##  [321,] 1.0525829
##  [322,] 1.0550626
##  [323,] 0.9589140
##  [324,] 1.0406578
##  [325,] 1.0322118
##  [326,] 0.9977062
##  [327,] 1.0039366
##  [328,] 0.9830758
##  [329,] 0.9921908
##  [330,] 0.9901266
##  [331,] 0.9071621
##  [332,] 1.0581722
##  [333,] 0.9275441
##  [334,] 0.9757797
##  [335,] 1.0364004
##  [336,] 1.0093758
##  [337,] 1.0925612
##  [338,] 0.9292176
##  [339,] 1.0447534
##  [340,] 1.0463701
##  [341,] 0.9663561
##  [342,] 1.1294943
##  [343,] 1.0543110
##  [344,] 0.9955550
##  [345,] 1.0290246
##  [346,] 1.0069820
##  [347,] 1.0133262
##  [348,] 1.0026611
##  [349,] 0.9743135
##  [350,] 0.9034668
##  [351,] 0.9787837
##  [352,] 0.9931086
##  [353,] 0.9495112
##  [354,] 1.0203716
##  [355,] 0.9633599
##  [356,] 0.9932725
##  [357,] 0.9188411
##  [358,] 0.9797156
##  [359,] 1.0263012
##  [360,] 1.0810994
##  [361,] 0.9122455
##  [362,] 0.9637441
##  [363,] 0.9975403
##  [364,] 0.9470584
##  [365,] 0.9641211
##  [366,] 0.9502657
##  [367,] 0.9304827
##  [368,] 1.1213904
##  [369,] 1.0195115
##  [370,] 0.9812633
##  [371,] 0.9945350
##  [372,] 1.0842029
##  [373,] 0.9800262
##  [374,] 1.0419747
##  [375,] 0.9593313
##  [376,] 1.0118931
##  [377,] 0.9636587
##  [378,] 0.9777803
##  [379,] 0.9844076
##  [380,] 1.0237296
##  [381,] 0.9676309
##  [382,] 0.9907608
##  [383,] 0.9509327
##  [384,] 1.0104312
##  [385,] 0.9898010
##  [386,] 0.9914596
##  [387,] 0.9926964
##  [388,] 0.9272390
##  [389,] 1.0859544
##  [390,] 0.9391759
##  [391,] 1.0646023
##  [392,] 0.9504514
##  [393,] 1.0404808
##  [394,] 0.9793928
##  [395,] 0.9630236
##  [396,] 0.9496162
##  [397,] 0.9764154
##  [398,] 1.0405933
##  [399,] 0.9598280
##  [400,] 0.9520389
##  [401,] 0.9312467
##  [402,] 1.0217220
##  [403,] 0.9904183
##  [404,] 1.0346964
##  [405,] 1.0100569
##  [406,] 0.9765250
##  [407,] 1.0579816
##  [408,] 0.9326834
##  [409,] 1.0585832
##  [410,] 0.9633531
##  [411,] 1.0005675
##  [412,] 0.9808573
##  [413,] 1.0300796
##  [414,] 0.9987657
##  [415,] 1.0037288
##  [416,] 0.9757207
##  [417,] 0.9646587
##  [418,] 1.0699999
##  [419,] 1.0165967
##  [420,] 1.0383626
##  [421,] 0.9703039
##  [422,] 0.9988245
##  [423,] 1.0962656
##  [424,] 0.9537702
##  [425,] 0.9892349
##  [426,] 0.9704986
##  [427,] 0.9613857
##  [428,] 1.0477329
##  [429,] 1.1016799
##  [430,] 0.9377356
##  [431,] 0.9044613
##  [432,] 1.0715584
##  [433,] 1.0838786
##  [434,] 0.9969844
##  [435,] 0.9326715
##  [436,] 0.9957885
##  [437,] 1.0065746
##  [438,] 1.1051500
##  [439,] 0.9709455
##  [440,] 1.0348427
##  [441,] 1.0012286
##  [442,] 0.9883767
##  [443,] 0.9261473
##  [444,] 1.0326658
##  [445,] 0.9287717
##  [446,] 0.9800595
##  [447,] 1.0553603
##  [448,] 1.0522890
##  [449,] 0.9757999
##  [450,] 1.0518046
##  [451,] 0.9921020
##  [452,] 0.9754764
##  [453,] 1.0298858
##  [454,] 0.9408529
##  [455,] 1.0250689
##  [456,] 1.0703684
##  [457,] 1.0408168
##  [458,] 0.9835806
##  [459,] 1.0629173
##  [460,] 1.0189364
##  [461,] 0.9780001
##  [462,] 1.0587656
##  [463,] 1.0044898
##  [464,] 1.0325159
##  [465,] 0.9089409
##  [466,] 0.9460213
##  [467,] 1.0118896
##  [468,] 0.9581067
##  [469,] 0.9960289
##  [470,] 0.9927821
##  [471,] 1.0375680
##  [472,] 1.1080328
##  [473,] 0.9595645
##  [474,] 0.9953810
##  [475,] 1.0260683
##  [476,] 0.9575311
##  [477,] 1.0216294
##  [478,] 0.9639473
##  [479,] 0.9864290
##  [480,] 0.9540851
##  [481,] 1.0804662
##  [482,] 1.0453476
##  [483,] 0.9948295
##  [484,] 1.0323705
##  [485,] 1.0023869
##  [486,] 0.9822555
##  [487,] 0.9565666
##  [488,] 0.9762558
##  [489,] 1.0575329
##  [490,] 1.0824875
##  [491,] 0.9865648
##  [492,] 0.9329112
##  [493,] 0.9900931
##  [494,] 0.9693573
##  [495,] 0.9636880
##  [496,] 0.9713056
##  [497,] 0.9901984
##  [498,] 1.0035762
##  [499,] 1.0995724
##  [500,] 0.8872206
##  [501,] 0.9910598
##  [502,] 0.9692349
##  [503,] 0.9325651
##  [504,] 1.0149416
##  [505,] 0.9763938
##  [506,] 1.0031036
##  [507,] 1.0081550
##  [508,] 1.0418998
##  [509,] 0.9995311
##  [510,] 1.0360653
##  [511,] 0.9762990
##  [512,] 1.0431633
##  [513,] 0.8996384
##  [514,] 0.9879828
##  [515,] 1.1147640
##  [516,] 0.9792805
##  [517,] 1.0105344
##  [518,] 1.0051141
##  [519,] 0.9798845
##  [520,] 0.9018919
##  [521,] 1.0295709
##  [522,] 0.9539894
##  [523,] 1.1044742
##  [524,] 1.0864593
##  [525,] 1.0524631
##  [526,] 0.9190339
##  [527,] 1.0230166
##  [528,] 1.0161344
##  [529,] 0.9716363
##  [530,] 1.0245420
##  [531,] 0.9438998
##  [532,] 0.9192135
##  [533,] 1.0207479
##  [534,] 0.9139306
##  [535,] 1.0271254
##  [536,] 1.0307403
##  [537,] 1.0047974
##  [538,] 1.1192621
##  [539,] 0.8946642
##  [540,] 1.1098444
##  [541,] 1.0092146
##  [542,] 0.9766015
##  [543,] 0.9965279
##  [544,] 0.9537294
##  [545,] 1.0923511
##  [546,] 0.9696849
##  [547,] 1.0356573
##  [548,] 0.9589194
##  [549,] 0.9669964
##  [550,] 0.9810144
##  [551,] 1.0017671
##  [552,] 0.9967778
##  [553,] 1.0041017
##  [554,] 1.0190893
##  [555,] 1.0262176
##  [556,] 1.0552424
##  [557,] 1.1373110
##  [558,] 1.0177691
##  [559,] 0.9185638
##  [560,] 0.9692808
##  [561,] 1.0119048
##  [562,] 0.9861838
##  [563,] 0.9528692
##  [564,] 1.0217292
##  [565,] 0.9617482
##  [566,] 1.0032463
##  [567,] 1.0421687
##  [568,] 1.1572612
##  [569,] 0.9646609
##  [570,] 0.9810458
##  [571,] 1.0277672
##  [572,] 0.9143013
##  [573,] 1.0054449
##  [574,] 1.0902208
##  [575,] 1.0101314
##  [576,] 1.0280450
##  [577,] 1.0373086
##  [578,] 0.9581777
##  [579,] 0.9335488
##  [580,] 0.9812232
##  [581,] 0.9883121
##  [582,] 0.9570692
##  [583,] 0.9214992
##  [584,] 0.9376947
##  [585,] 1.0234933
##  [586,] 0.9987106
##  [587,] 0.9619668
##  [588,] 1.0421421
##  [589,] 0.9616677
##  [590,] 1.0424010
##  [591,] 1.0317432
##  [592,] 0.9848182
##  [593,] 1.0498072
##  [594,] 0.9546565
##  [595,] 0.9728978
##  [596,] 0.8978994
##  [597,] 0.9993725
##  [598,] 0.9980583
##  [599,] 0.9144635
##  [600,] 1.0033382
##  [601,] 0.9630091
##  [602,] 1.0852917
##  [603,] 0.9470832
##  [604,] 1.0983455
##  [605,] 1.0602178
##  [606,] 1.1002114
##  [607,] 0.9109027
##  [608,] 1.0647928
##  [609,] 0.9523116
##  [610,] 0.9261020
##  [611,] 0.9538194
##  [612,] 0.9543791
##  [613,] 1.0458819
##  [614,] 0.9823375
##  [615,] 1.0142394
##  [616,] 0.9065677
##  [617,] 1.0614596
##  [618,] 1.0225937
##  [619,] 1.1314475
##  [620,] 1.0930800
##  [621,] 1.0986326
##  [622,] 1.0120093
##  [623,] 1.0315417
##  [624,] 1.0142423
##  [625,] 0.9036410
##  [626,] 0.9803686
##  [627,] 1.0015576
##  [628,] 0.9654637
##  [629,] 1.0023089
##  [630,] 0.9107155
##  [631,] 0.9503778
##  [632,] 1.0551418
##  [633,] 0.9165015
##  [634,] 0.9099992
##  [635,] 0.9931081
##  [636,] 1.0213229
##  [637,] 1.0497030
##  [638,] 0.9576490
##  [639,] 1.0218758
##  [640,] 1.0684583
##  [641,] 0.9128959
##  [642,] 1.0098816
##  [643,] 0.9695580
##  [644,] 0.9122037
##  [645,] 1.0210078
##  [646,] 0.8777000
##  [647,] 0.9886205
##  [648,] 0.9691386
##  [649,] 0.9974739
##  [650,] 0.9389521
##  [651,] 1.0031312
##  [652,] 0.8944307
##  [653,] 0.9163596
##  [654,] 0.9821764
##  [655,] 0.9931933
##  [656,] 1.0017546
##  [657,] 1.0697733
##  [658,] 1.0331394
##  [659,] 1.0114286
##  [660,] 1.0840469
##  [661,] 1.0019993
##  [662,] 0.9965815
##  [663,] 1.1058227
##  [664,] 1.0327332
##  [665,] 0.9298871
##  [666,] 1.0426173
##  [667,] 0.9954298
##  [668,] 1.0029417
##  [669,] 0.9529008
##  [670,] 1.0098017
##  [671,] 1.0029503
##  [672,] 1.0895663
##  [673,] 0.9905838
##  [674,] 0.9833297
##  [675,] 1.0472262
##  [676,] 1.0794736
##  [677,] 0.9986070
##  [678,] 0.9584243
##  [679,] 0.9653499
##  [680,] 1.0477719
##  [681,] 1.0010006
##  [682,] 0.9135224
##  [683,] 0.9927375
##  [684,] 0.9630487
##  [685,] 0.9846343
##  [686,] 1.0067561
##  [687,] 1.1145442
##  [688,] 0.9777353
##  [689,] 1.0906787
##  [690,] 0.9319301
##  [691,] 0.9710592
##  [692,] 0.9660963
##  [693,] 1.0592910
##  [694,] 1.0555848
##  [695,] 1.0276567
##  [696,] 0.9468888
##  [697,] 1.0859180
##  [698,] 1.0098352
##  [699,] 1.0195330
##  [700,] 1.0164622
##  [701,] 1.0188161
##  [702,] 1.0228751
##  [703,] 0.9867352
##  [704,] 1.0701812
##  [705,] 1.0399874
##  [706,] 0.9483412
##  [707,] 1.0270869
##  [708,] 0.9639333
##  [709,] 1.0404621
##  [710,] 1.1285062
##  [711,] 1.0265614
##  [712,] 0.9224290
##  [713,] 0.9679860
##  [714,] 0.9636457
##  [715,] 0.9897742
##  [716,] 0.9453435
##  [717,] 1.0258043
##  [718,] 1.1215423
##  [719,] 0.9334081
##  [720,] 0.9883168
##  [721,] 0.9596312
##  [722,] 0.9853092
##  [723,] 1.0967124
##  [724,] 0.9193019
##  [725,] 1.0675125
##  [726,] 0.9914393
##  [727,] 0.9578349
##  [728,] 1.0043166
##  [729,] 1.0602241
##  [730,] 1.0336570
##  [731,] 0.9850322
##  [732,] 0.9734680
##  [733,] 1.0467124
##  [734,] 0.9148257
##  [735,] 0.9353717
##  [736,] 0.9273836
##  [737,] 0.9773583
##  [738,] 0.9881600
##  [739,] 0.9540619
##  [740,] 0.9748602
##  [741,] 0.9450478
##  [742,] 1.0281136
##  [743,] 0.9606092
##  [744,] 1.0431225
##  [745,] 1.0843065
##  [746,] 0.9694639
##  [747,] 1.1049766
##  [748,] 1.0084460
##  [749,] 1.0564636
##  [750,] 1.0365434
##  [751,] 0.9973700
##  [752,] 0.9453330
##  [753,] 1.0369624
##  [754,] 1.0219429
##  [755,] 0.8928597
##  [756,] 0.9820080
##  [757,] 1.0854704
##  [758,] 1.0722857
##  [759,] 0.9356478
##  [760,] 0.9896495
##  [761,] 1.0303169
##  [762,] 0.9275621
##  [763,] 1.0086027
##  [764,] 0.9177014
##  [765,] 0.9634274
##  [766,] 1.0248965
##  [767,] 0.9734796
##  [768,] 1.0155049
##  [769,] 1.0030661
##  [770,] 0.8998902
##  [771,] 1.0334104
##  [772,] 1.0365089
##  [773,] 0.9854934
##  [774,] 0.9337654
##  [775,] 1.0085557
##  [776,] 0.9679819
##  [777,] 0.9196501
##  [778,] 1.0485212
##  [779,] 1.0058944
##  [780,] 0.9793906
##  [781,] 1.0122181
##  [782,] 0.9885109
##  [783,] 1.0020738
##  [784,] 1.0084445
##  [785,] 1.0343408
##  [786,] 0.9374885
##  [787,] 1.0710238
##  [788,] 0.9568121
##  [789,] 1.0362300
##  [790,] 0.9829519
##  [791,] 0.9805421
##  [792,] 0.9709190
##  [793,] 1.0332069
##  [794,] 0.9824008
##  [795,] 0.9522221
##  [796,] 0.9495027
##  [797,] 0.9581032
##  [798,] 1.0117621
##  [799,] 0.9942186
##  [800,] 1.0099426
##  [801,] 0.9536613
##  [802,] 0.9507380
##  [803,] 0.9632088
##  [804,] 1.0193233
##  [805,] 0.9975944
##  [806,] 1.0292642
##  [807,] 0.9363680
##  [808,] 0.9490121
##  [809,] 0.9321500
##  [810,] 0.9921772
##  [811,] 0.9221894
##  [812,] 1.0139916
##  [813,] 1.0076909
##  [814,] 0.9237817
##  [815,] 1.0218752
##  [816,] 1.0432873
##  [817,] 0.9967735
##  [818,] 0.9752456
##  [819,] 0.9694460
##  [820,] 1.0820693
##  [821,] 1.0353072
##  [822,] 1.0062434
##  [823,] 0.9669808
##  [824,] 0.9961568
##  [825,] 1.0851365
##  [826,] 1.1511599
##  [827,] 1.0220628
##  [828,] 1.0561518
##  [829,] 0.9629563
##  [830,] 1.0359038
##  [831,] 0.9792914
##  [832,] 0.9371576
##  [833,] 0.9778148
##  [834,] 1.0037138
##  [835,] 1.1364375
##  [836,] 0.9910399
##  [837,] 0.9265723
##  [838,] 0.9879833
##  [839,] 1.0755498
##  [840,] 0.8835080
##  [841,] 1.0856037
##  [842,] 0.9014485
##  [843,] 0.9948904
##  [844,] 1.0026498
##  [845,] 0.9715380
##  [846,] 0.9759157
##  [847,] 0.9449908
##  [848,] 1.0143413
##  [849,] 0.9629399
##  [850,] 1.0818777
##  [851,] 1.0888577
##  [852,] 1.0388463
##  [853,] 1.0100195
##  [854,] 0.9911224
##  [855,] 1.0343968
##  [856,] 0.9913102
##  [857,] 1.0471214
##  [858,] 0.9911623
##  [859,] 1.1036716
##  [860,] 1.0102292
##  [861,] 1.0481753
##  [862,] 0.9920415
##  [863,] 1.0612221
##  [864,] 1.1317762
##  [865,] 0.9434227
##  [866,] 0.9762293
##  [867,] 1.0185329
##  [868,] 1.0392037
##  [869,] 1.0303036
##  [870,] 0.9921429
##  [871,] 0.9911166
##  [872,] 0.8849189
##  [873,] 0.9468581
##  [874,] 0.9516668
##  [875,] 1.0206849
##  [876,] 0.9690388
##  [877,] 0.9843756
##  [878,] 1.1167040
##  [879,] 0.9918838
##  [880,] 0.9796883
##  [881,] 1.0199871
##  [882,] 0.9944763
##  [883,] 1.0300061
##  [884,] 1.0473769
##  [885,] 1.0577851
##  [886,] 1.1133129
##  [887,] 0.9406189
##  [888,] 0.8860257
##  [889,] 1.0718064
##  [890,] 1.0100385
##  [891,] 1.0838592
##  [892,] 1.0823296
##  [893,] 1.1110922
##  [894,] 0.9987832
##  [895,] 0.9998272
##  [896,] 1.0031388
##  [897,] 1.0736002
##  [898,] 1.1128228
##  [899,] 0.9038780
##  [900,] 1.1087672
##  [901,] 0.9788629
##  [902,] 1.0436966
##  [903,] 0.8900487
##  [904,] 0.9319297
##  [905,] 0.9886329
##  [906,] 1.0628822
##  [907,] 0.9960263
##  [908,] 0.9795310
##  [909,] 0.9794464
##  [910,] 1.0350590
##  [911,] 1.1215992
##  [912,] 1.0080596
##  [913,] 1.0001647
##  [914,] 1.0342052
##  [915,] 0.9901320
##  [916,] 0.9827556
##  [917,] 0.8986972
##  [918,] 0.9763659
##  [919,] 1.0237301
##  [920,] 1.0783971
##  [921,] 0.9860015
##  [922,] 1.0601323
##  [923,] 0.9264739
##  [924,] 1.0813827
##  [925,] 1.0234906
##  [926,] 0.9700150
##  [927,] 1.0157215
##  [928,] 1.0099746
##  [929,] 1.0256146
##  [930,] 1.0090254
##  [931,] 0.9539743
##  [932,] 1.0395800
##  [933,] 1.0744417
##  [934,] 0.9302110
##  [935,] 0.9729301
##  [936,] 1.0103567
##  [937,] 1.0194747
##  [938,] 1.0053909
##  [939,] 1.0266038
##  [940,] 1.0689224
##  [941,] 1.0361066
##  [942,] 1.0386093
##  [943,] 1.0907035
##  [944,] 1.0265465
##  [945,] 0.9980115
##  [946,] 1.0200651
##  [947,] 1.0247600
##  [948,] 1.1410642
##  [949,] 0.9479641
##  [950,] 0.9832687
##  [951,] 1.0018967
##  [952,] 1.1147367
##  [953,] 0.9416490
##  [954,] 1.0160332
##  [955,] 0.9460230
##  [956,] 0.9515729
##  [957,] 1.0712917
##  [958,] 1.0384385
##  [959,] 0.9296148
##  [960,] 0.9887829
##  [961,] 1.0140426
##  [962,] 0.9700145
##  [963,] 1.0501663
##  [964,] 1.0362480
##  [965,] 0.9915232
##  [966,] 1.0677256
##  [967,] 1.0252473
##  [968,] 1.0376591
##  [969,] 1.1392649
##  [970,] 1.0364227
##  [971,] 0.9386362
##  [972,] 1.0759466
##  [973,] 0.9390613
##  [974,] 1.0467323
##  [975,] 0.9869115
##  [976,] 1.1346822
##  [977,] 1.0155135
##  [978,] 1.0260794
##  [979,] 0.9616699
##  [980,] 1.0429225
##  [981,] 1.0016666
##  [982,] 0.9954694
##  [983,] 1.0477096
##  [984,] 0.9608299
##  [985,] 0.9942137
##  [986,] 1.0286936
##  [987,] 1.0119473
##  [988,] 1.0165418
##  [989,] 0.9571785
##  [990,] 1.0046239
##  [991,] 0.9413403
##  [992,] 0.9402256
##  [993,] 0.9160779
##  [994,] 0.9155975
##  [995,] 0.9995359
##  [996,] 1.0472588
##  [997,] 0.9866608
##  [998,] 0.9680571
##  [999,] 0.9488023
## 
## $model.matrix
##    (Intercept) site_pH1 site_pH2 site_pH3 site_pH4
## 1            1        1        0        0        0
## 2            1        0        1        0        0
## 3            1        0        0        1        0
## 4            1        0        0        0        1
## 5            1       -1       -1       -1       -1
## 6            1        1        0        0        0
## 7            1        0        1        0        0
## 8            1        0        0        1        0
## 9            1        0        0        0        1
## 10           1       -1       -1       -1       -1
## 11           1        1        0        0        0
## 12           1        0        1        0        0
## 13           1        0        0        1        0
## 14           1        0        0        0        1
## 15           1       -1       -1       -1       -1
## 16           1        1        0        0        0
## 17           1        0        1        0        0
## 18           1        0        0        1        0
## 19           1        0        0        0        1
## 20           1       -1       -1       -1       -1
## 
## $terms
## dist_tab_assay ~ site_pH
## attr(,"variables")
## list(dist_tab_assay, site_pH)
## attr(,"factors")
##                site_pH
## dist_tab_assay       0
## site_pH              1
## attr(,"term.labels")
## [1] "site_pH"
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
```

```r
anova(betadisper(dist_tab_assay,samples3$site_pH))
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df  Sum Sq Mean Sq F value Pr(>F)
## Groups     4  7121.9  1780.5  1.1231 0.3825
## Residuals 15 23780.2  1585.3
```

```r
count_tab_assay <- assay(vsdNatSim)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samplesNatSim,dist_tab_assay ~ site, method="euclidian")
```

```
## 'adonis' will be deprecated: use 'adonis2' instead
```

```
## $aov.tab
## Permutation: free
## Number of permutations: 999
## 
## Terms added sequentially (first to last)
## 
##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
## site       1     35434   35434  1.1619 0.16223  0.028 *
## Residuals  6    182977   30496         0.83777         
## Total      7    218411                 1.00000         
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## $call
## adonis(formula = dist_tab_assay ~ site, data = samplesNatSim, 
##     method = "euclidian")
## 
## $coefficients
## NULL
## 
## $coef.sites
##                   [,1]      [,2]      [,3]      [,4]       [,5]      [,6]       [,7]      [,8]
## (Intercept) 182.458863 180.14862 351.77699 176.43422 181.113773 174.36441 178.799745 176.63663
## site1         2.743034  56.59463 -49.84527  52.64912   6.165957  56.00634   6.800219  55.45145
## 
## $f.perms
##              [,1]
##    [1,] 1.0523772
##    [2,] 1.1618997
##    [3,] 0.9832302
##    [4,] 0.9674016
##    [5,] 0.9319773
##    [6,] 1.0844318
##    [7,] 1.0427419
##    [8,] 0.9049066
##    [9,] 1.0645854
##   [10,] 1.0427419
##   [11,] 0.9976266
##   [12,] 0.9237132
##   [13,] 1.1618997
##   [14,] 0.9542853
##   [15,] 0.9185880
##   [16,] 1.0273393
##   [17,] 1.0000467
##   [18,] 1.1056910
##   [19,] 1.0645854
##   [20,] 1.0000467
##   [21,] 1.0523772
##   [22,] 1.0085864
##   [23,] 0.9965271
##   [24,] 0.9319773
##   [25,] 1.0028665
##   [26,] 0.9976266
##   [27,] 1.0273393
##   [28,] 0.9858118
##   [29,] 0.9361535
##   [30,] 0.9965271
##   [31,] 1.0273393
##   [32,] 0.9237132
##   [33,] 0.9674016
##   [34,] 0.9828184
##   [35,] 0.9858118
##   [36,] 1.0085864
##   [37,] 0.9185880
##   [38,] 0.9237132
##   [39,] 0.9944049
##   [40,] 0.9185880
##   [41,] 0.9185880
##   [42,] 1.1618997
##   [43,] 0.9832302
##   [44,] 0.9319773
##   [45,] 1.0039037
##   [46,] 0.9681868
##   [47,] 0.9746116
##   [48,] 1.1618997
##   [49,] 0.9049066
##   [50,] 0.9049066
##   [51,] 0.9319773
##   [52,] 0.9944049
##   [53,] 0.9681868
##   [54,] 1.0844318
##   [55,] 0.9144307
##   [56,] 0.9049066
##   [57,] 1.0028665
##   [58,] 1.1056910
##   [59,] 1.0523772
##   [60,] 0.9144307
##   [61,] 0.9944049
##   [62,] 1.1056910
##   [63,] 1.0523772
##   [64,] 0.9462140
##   [65,] 0.9319773
##   [66,] 0.9049066
##   [67,] 1.1006436
##   [68,] 0.9319773
##   [69,] 1.1006436
##   [70,] 0.9828184
##   [71,] 1.0523772
##   [72,] 0.9944049
##   [73,] 0.9944049
##   [74,] 0.9361535
##   [75,] 0.9185880
##   [76,] 0.9144307
##   [77,] 0.9681868
##   [78,] 0.9319773
##   [79,] 1.0028665
##   [80,] 1.0085864
##   [81,] 0.9944049
##   [82,] 0.9965271
##   [83,] 1.0085864
##   [84,] 1.0039037
##   [85,] 0.9828184
##   [86,] 0.9746116
##   [87,] 1.1618997
##   [88,] 0.9361535
##   [89,] 1.0758820
##   [90,] 1.0427419
##   [91,] 1.1618997
##   [92,] 0.9542853
##   [93,] 0.9185880
##   [94,] 1.0575953
##   [95,] 1.0085864
##   [96,] 0.9808044
##   [97,] 1.0758820
##   [98,] 0.9965271
##   [99,] 0.9681868
##  [100,] 1.0655676
##  [101,] 0.9319773
##  [102,] 0.9858118
##  [103,] 1.0039037
##  [104,] 0.9674016
##  [105,] 0.9319773
##  [106,] 0.9976266
##  [107,] 1.0575953
##  [108,] 0.9746116
##  [109,] 0.9319773
##  [110,] 0.9319773
##  [111,] 1.0427419
##  [112,] 1.0844318
##  [113,] 0.9828184
##  [114,] 0.9185880
##  [115,] 0.9361535
##  [116,] 1.0575953
##  [117,] 0.9237132
##  [118,] 0.9134597
##  [119,] 0.9746116
##  [120,] 0.9185880
##  [121,] 1.0758820
##  [122,] 0.9976266
##  [123,] 1.0523772
##  [124,] 1.0039037
##  [125,] 0.9746116
##  [126,] 0.9361535
##  [127,] 0.9828184
##  [128,] 0.9808044
##  [129,] 0.9674016
##  [130,] 1.0085864
##  [131,] 0.9808044
##  [132,] 0.9237132
##  [133,] 1.0645854
##  [134,] 1.0575953
##  [135,] 1.0039037
##  [136,] 0.9746116
##  [137,] 0.9144307
##  [138,] 0.9319773
##  [139,] 0.9832302
##  [140,] 0.9894063
##  [141,] 0.9185880
##  [142,] 1.0000467
##  [143,] 1.0085864
##  [144,] 1.0575953
##  [145,] 0.9319773
##  [146,] 1.0000467
##  [147,] 1.0085864
##  [148,] 1.0085864
##  [149,] 0.9681868
##  [150,] 0.9965271
##  [151,] 0.9828184
##  [152,] 0.9462140
##  [153,] 0.9976266
##  [154,] 1.1006436
##  [155,] 1.0028665
##  [156,] 1.0273393
##  [157,] 1.0523772
##  [158,] 0.9965271
##  [159,] 1.0575953
##  [160,] 1.0039037
##  [161,] 1.0085864
##  [162,] 0.9134597
##  [163,] 0.9944049
##  [164,] 1.1618997
##  [165,] 0.9681868
##  [166,] 0.9462140
##  [167,] 0.9361535
##  [168,] 0.9049066
##  [169,] 1.1006436
##  [170,] 1.1618997
##  [171,] 0.9965271
##  [172,] 0.9858118
##  [173,] 0.9858118
##  [174,] 0.9674016
##  [175,] 1.0844318
##  [176,] 0.9894063
##  [177,] 0.9361535
##  [178,] 1.1618997
##  [179,] 1.0655676
##  [180,] 0.9319773
##  [181,] 0.9894063
##  [182,] 0.9858118
##  [183,] 0.9681868
##  [184,] 0.9462140
##  [185,] 0.9828184
##  [186,] 0.9462140
##  [187,] 0.9237132
##  [188,] 1.0645854
##  [189,] 1.0000467
##  [190,] 1.0758820
##  [191,] 0.9144307
##  [192,] 0.9976266
##  [193,] 0.9319773
##  [194,] 0.9237132
##  [195,] 1.0758820
##  [196,] 0.9542853
##  [197,] 0.9894063
##  [198,] 0.9319773
##  [199,] 0.9894063
##  [200,] 0.9319773
##  [201,] 0.9185880
##  [202,] 1.0844318
##  [203,] 1.0028665
##  [204,] 0.9858118
##  [205,] 0.9319773
##  [206,] 1.0655676
##  [207,] 1.0085864
##  [208,] 1.1056910
##  [209,] 0.9134597
##  [210,] 0.9049066
##  [211,] 0.9049066
##  [212,] 0.9832302
##  [213,] 0.9832302
##  [214,] 0.9049066
##  [215,] 0.9746116
##  [216,] 1.0028665
##  [217,] 1.1618997
##  [218,] 0.9894063
##  [219,] 1.0427419
##  [220,] 1.0000467
##  [221,] 1.0758820
##  [222,] 0.9319773
##  [223,] 1.0039037
##  [224,] 0.9894063
##  [225,] 1.0575953
##  [226,] 0.9858118
##  [227,] 1.0028665
##  [228,] 1.0000467
##  [229,] 1.0028665
##  [230,] 0.9462140
##  [231,] 1.0655676
##  [232,] 1.0000467
##  [233,] 0.9542853
##  [234,] 0.9237132
##  [235,] 0.9828184
##  [236,] 0.9828184
##  [237,] 1.0523772
##  [238,] 0.9049066
##  [239,] 1.1006436
##  [240,] 0.9681868
##  [241,] 0.9746116
##  [242,] 0.9965271
##  [243,] 1.1056910
##  [244,] 1.1618997
##  [245,] 1.0427419
##  [246,] 1.0085864
##  [247,] 1.0645854
##  [248,] 1.0523772
##  [249,] 1.0645854
##  [250,] 1.0273393
##  [251,] 1.0575953
##  [252,] 0.9049066
##  [253,] 0.9144307
##  [254,] 0.9858118
##  [255,] 0.9542853
##  [256,] 1.0427419
##  [257,] 1.0039037
##  [258,] 0.9808044
##  [259,] 1.1056910
##  [260,] 0.9681868
##  [261,] 1.0028665
##  [262,] 0.9144307
##  [263,] 0.9462140
##  [264,] 1.0645854
##  [265,] 0.9965271
##  [266,] 1.0523772
##  [267,] 1.0085864
##  [268,] 0.9976266
##  [269,] 1.0655676
##  [270,] 1.0844318
##  [271,] 1.0000467
##  [272,] 1.0039037
##  [273,] 1.0427419
##  [274,] 0.9976266
##  [275,] 1.0028665
##  [276,] 0.9319773
##  [277,] 1.0523772
##  [278,] 1.0273393
##  [279,] 1.0000467
##  [280,] 1.0039037
##  [281,] 1.0039037
##  [282,] 0.9828184
##  [283,] 0.9976266
##  [284,] 0.9462140
##  [285,] 0.9674016
##  [286,] 0.9681868
##  [287,] 1.0844318
##  [288,] 0.9746116
##  [289,] 0.9319773
##  [290,] 0.9144307
##  [291,] 0.9746116
##  [292,] 0.9185880
##  [293,] 0.9319773
##  [294,] 0.9674016
##  [295,] 1.0523772
##  [296,] 0.9462140
##  [297,] 1.0427419
##  [298,] 0.9858118
##  [299,] 1.0039037
##  [300,] 0.9894063
##  [301,] 0.9681868
##  [302,] 1.1006436
##  [303,] 0.9361535
##  [304,] 0.9049066
##  [305,] 0.9144307
##  [306,] 1.0028665
##  [307,] 0.9185880
##  [308,] 0.9894063
##  [309,] 0.9049066
##  [310,] 1.0645854
##  [311,] 0.9542853
##  [312,] 1.0523772
##  [313,] 1.0575953
##  [314,] 0.9361535
##  [315,] 1.1006436
##  [316,] 1.0000467
##  [317,] 1.0758820
##  [318,] 1.0273393
##  [319,] 0.9542853
##  [320,] 0.9542853
##  [321,] 1.0523772
##  [322,] 1.0085864
##  [323,] 0.9976266
##  [324,] 0.9144307
##  [325,] 0.9542853
##  [326,] 0.9049066
##  [327,] 0.9832302
##  [328,] 1.0523772
##  [329,] 1.0523772
##  [330,] 0.9134597
##  [331,] 0.9828184
##  [332,] 1.0523772
##  [333,] 1.1056910
##  [334,] 0.9894063
##  [335,] 0.9808044
##  [336,] 0.9049066
##  [337,] 1.0273393
##  [338,] 1.0844318
##  [339,] 0.9319773
##  [340,] 1.0000467
##  [341,] 0.9965271
##  [342,] 1.0039037
##  [343,] 1.0645854
##  [344,] 0.9462140
##  [345,] 0.9319773
##  [346,] 1.0575953
##  [347,] 0.9828184
##  [348,] 0.9976266
##  [349,] 0.9894063
##  [350,] 1.0523772
##  [351,] 1.0085864
##  [352,] 0.9542853
##  [353,] 1.0844318
##  [354,] 1.0758820
##  [355,] 1.1006436
##  [356,] 1.0844318
##  [357,] 1.1056910
##  [358,] 0.9185880
##  [359,] 1.0758820
##  [360,] 1.0844318
##  [361,] 1.0645854
##  [362,] 0.9808044
##  [363,] 0.9828184
##  [364,] 0.9237132
##  [365,] 0.9944049
##  [366,] 1.0039037
##  [367,] 0.9965271
##  [368,] 0.9462140
##  [369,] 1.0039037
##  [370,] 1.1618997
##  [371,] 1.0427419
##  [372,] 0.9144307
##  [373,] 0.9674016
##  [374,] 1.0575953
##  [375,] 0.9185880
##  [376,] 1.0655676
##  [377,] 1.1056910
##  [378,] 1.0844318
##  [379,] 1.0655676
##  [380,] 1.0427419
##  [381,] 0.9808044
##  [382,] 1.0028665
##  [383,] 0.9858118
##  [384,] 0.9808044
##  [385,] 0.9361535
##  [386,] 0.9542853
##  [387,] 0.9965271
##  [388,] 0.9674016
##  [389,] 0.9965271
##  [390,] 0.9462140
##  [391,] 1.0427419
##  [392,] 1.0028665
##  [393,] 0.9681868
##  [394,] 0.9542853
##  [395,] 0.9858118
##  [396,] 1.0844318
##  [397,] 1.0523772
##  [398,] 1.1006436
##  [399,] 1.0427419
##  [400,] 1.0575953
##  [401,] 1.0028665
##  [402,] 0.9674016
##  [403,] 0.9674016
##  [404,] 1.1056910
##  [405,] 1.0427419
##  [406,] 0.9462140
##  [407,] 1.0575953
##  [408,] 0.9462140
##  [409,] 0.9858118
##  [410,] 0.9462140
##  [411,] 1.0655676
##  [412,] 0.9237132
##  [413,] 1.0039037
##  [414,] 1.0000467
##  [415,] 0.9965271
##  [416,] 1.0039037
##  [417,] 0.9237132
##  [418,] 0.9944049
##  [419,] 1.0655676
##  [420,] 1.0028665
##  [421,] 0.9894063
##  [422,] 1.0645854
##  [423,] 0.9976266
##  [424,] 0.9828184
##  [425,] 0.9462140
##  [426,] 1.0000467
##  [427,] 1.0039037
##  [428,] 0.9746116
##  [429,] 1.0039037
##  [430,] 0.9144307
##  [431,] 1.0523772
##  [432,] 1.0273393
##  [433,] 0.9808044
##  [434,] 0.9462140
##  [435,] 0.9049066
##  [436,] 1.0645854
##  [437,] 0.9049066
##  [438,] 1.1618997
##  [439,] 0.9894063
##  [440,] 1.0844318
##  [441,] 0.9976266
##  [442,] 0.9542853
##  [443,] 0.9237132
##  [444,] 0.9894063
##  [445,] 1.0575953
##  [446,] 1.0523772
##  [447,] 0.9828184
##  [448,] 0.9828184
##  [449,] 1.1056910
##  [450,] 0.9462140
##  [451,] 1.0028665
##  [452,] 1.0575953
##  [453,] 0.9319773
##  [454,] 0.9674016
##  [455,] 1.0645854
##  [456,] 1.0085864
##  [457,] 1.1006436
##  [458,] 0.9185880
##  [459,] 0.9976266
##  [460,] 1.0085864
##  [461,] 0.9746116
##  [462,] 0.9185880
##  [463,] 0.9944049
##  [464,] 1.1618997
##  [465,] 1.0523772
##  [466,] 1.0273393
##  [467,] 1.0273393
##  [468,] 1.0645854
##  [469,] 0.9237132
##  [470,] 1.0000467
##  [471,] 1.0039037
##  [472,] 0.9319773
##  [473,] 1.1056910
##  [474,] 1.0085864
##  [475,] 1.0273393
##  [476,] 1.0273393
##  [477,] 0.9134597
##  [478,] 1.0758820
##  [479,] 1.0523772
##  [480,] 1.0427419
##  [481,] 0.9237132
##  [482,] 1.0000467
##  [483,] 1.0028665
##  [484,] 0.9832302
##  [485,] 0.9237132
##  [486,] 0.9361535
##  [487,] 1.0028665
##  [488,] 0.9319773
##  [489,] 1.0085864
##  [490,] 0.9049066
##  [491,] 0.9237132
##  [492,] 0.9185880
##  [493,] 1.1056910
##  [494,] 1.0645854
##  [495,] 1.0039037
##  [496,] 1.0575953
##  [497,] 1.0655676
##  [498,] 0.9049066
##  [499,] 0.9319773
##  [500,] 0.9134597
##  [501,] 0.9828184
##  [502,] 1.0028665
##  [503,] 0.9185880
##  [504,] 0.9858118
##  [505,] 0.9944049
##  [506,] 1.0758820
##  [507,] 1.0645854
##  [508,] 0.9894063
##  [509,] 1.0028665
##  [510,] 1.0000467
##  [511,] 0.9681868
##  [512,] 0.9681868
##  [513,] 0.9681868
##  [514,] 1.0273393
##  [515,] 0.9185880
##  [516,] 1.0085864
##  [517,] 0.9361535
##  [518,] 1.0844318
##  [519,] 1.1618997
##  [520,] 0.9746116
##  [521,] 1.0758820
##  [522,] 0.9542853
##  [523,] 0.9832302
##  [524,] 0.9361535
##  [525,] 0.9976266
##  [526,] 0.9976266
##  [527,] 0.9144307
##  [528,] 0.9185880
##  [529,] 1.0645854
##  [530,] 1.0028665
##  [531,] 0.9542853
##  [532,] 1.0000467
##  [533,] 1.0028665
##  [534,] 0.9237132
##  [535,] 0.9944049
##  [536,] 1.1006436
##  [537,] 1.0844318
##  [538,] 1.0000467
##  [539,] 0.9944049
##  [540,] 1.0039037
##  [541,] 0.9237132
##  [542,] 1.0844318
##  [543,] 0.9976266
##  [544,] 1.0085864
##  [545,] 0.9049066
##  [546,] 1.0758820
##  [547,] 1.0523772
##  [548,] 0.9542853
##  [549,] 0.9237132
##  [550,] 0.9894063
##  [551,] 1.0645854
##  [552,] 0.9894063
##  [553,] 0.9828184
##  [554,] 0.9976266
##  [555,] 0.9319773
##  [556,] 0.9976266
##  [557,] 0.9944049
##  [558,] 1.0523772
##  [559,] 1.0758820
##  [560,] 1.1056910
##  [561,] 0.9361535
##  [562,] 0.9976266
##  [563,] 0.9894063
##  [564,] 1.0655676
##  [565,] 0.9828184
##  [566,] 0.9237132
##  [567,] 1.1056910
##  [568,] 0.9681868
##  [569,] 1.0085864
##  [570,] 1.0655676
##  [571,] 1.0028665
##  [572,] 1.0645854
##  [573,] 0.9185880
##  [574,] 0.9237132
##  [575,] 1.1006436
##  [576,] 0.9858118
##  [577,] 0.9832302
##  [578,] 1.0273393
##  [579,] 0.9361535
##  [580,] 0.9237132
##  [581,] 0.9462140
##  [582,] 1.0655676
##  [583,] 0.9808044
##  [584,] 0.9746116
##  [585,] 0.9361535
##  [586,] 0.9049066
##  [587,] 1.0655676
##  [588,] 0.9674016
##  [589,] 1.1618997
##  [590,] 1.0039037
##  [591,] 1.0039037
##  [592,] 0.9976266
##  [593,] 0.9237132
##  [594,] 1.0427419
##  [595,] 0.9185880
##  [596,] 0.9832302
##  [597,] 1.0655676
##  [598,] 1.1618997
##  [599,] 1.0844318
##  [600,] 1.0575953
##  [601,] 0.9049066
##  [602,] 0.9185880
##  [603,] 1.0273393
##  [604,] 0.9185880
##  [605,] 0.9049066
##  [606,] 0.9237132
##  [607,] 1.0655676
##  [608,] 0.9858118
##  [609,] 0.9858118
##  [610,] 0.9681868
##  [611,] 0.9894063
##  [612,] 0.9808044
##  [613,] 1.0844318
##  [614,] 0.9858118
##  [615,] 1.0523772
##  [616,] 0.9746116
##  [617,] 1.0273393
##  [618,] 1.1618997
##  [619,] 1.0028665
##  [620,] 0.9542853
##  [621,] 0.9144307
##  [622,] 0.9746116
##  [623,] 0.9681868
##  [624,] 0.9965271
##  [625,] 0.9462140
##  [626,] 1.0085864
##  [627,] 1.1056910
##  [628,] 1.0039037
##  [629,] 1.0575953
##  [630,] 0.9361535
##  [631,] 1.0844318
##  [632,] 0.9049066
##  [633,] 1.0758820
##  [634,] 0.9828184
##  [635,] 1.1618997
##  [636,] 1.1056910
##  [637,] 0.9361535
##  [638,] 0.9134597
##  [639,] 1.0273393
##  [640,] 0.9944049
##  [641,] 0.9134597
##  [642,] 0.9976266
##  [643,] 0.9237132
##  [644,] 0.9049066
##  [645,] 0.9828184
##  [646,] 1.0758820
##  [647,] 0.9361535
##  [648,] 0.9944049
##  [649,] 1.0655676
##  [650,] 0.9144307
##  [651,] 0.9674016
##  [652,] 0.9237132
##  [653,] 1.0575953
##  [654,] 0.9237132
##  [655,] 0.9361535
##  [656,] 1.0427419
##  [657,] 0.9965271
##  [658,] 0.9976266
##  [659,] 0.9144307
##  [660,] 1.1006436
##  [661,] 1.0273393
##  [662,] 1.0427419
##  [663,] 0.9185880
##  [664,] 1.0655676
##  [665,] 0.9976266
##  [666,] 1.0028665
##  [667,] 0.9674016
##  [668,] 0.9049066
##  [669,] 1.1056910
##  [670,] 0.9894063
##  [671,] 1.0085864
##  [672,] 0.9542853
##  [673,] 0.9828184
##  [674,] 0.9828184
##  [675,] 0.9319773
##  [676,] 0.9462140
##  [677,] 1.0039037
##  [678,] 0.9832302
##  [679,] 1.0523772
##  [680,] 0.9832302
##  [681,] 0.9462140
##  [682,] 0.9681868
##  [683,] 0.9237132
##  [684,] 0.9976266
##  [685,] 0.9858118
##  [686,] 1.0758820
##  [687,] 1.0844318
##  [688,] 0.9185880
##  [689,] 1.0039037
##  [690,] 0.9976266
##  [691,] 1.0655676
##  [692,] 0.9144307
##  [693,] 0.9944049
##  [694,] 0.9858118
##  [695,] 1.0523772
##  [696,] 0.9681868
##  [697,] 0.9134597
##  [698,] 1.0273393
##  [699,] 1.0028665
##  [700,] 0.9944049
##  [701,] 0.9237132
##  [702,] 0.9542853
##  [703,] 1.0000467
##  [704,] 1.1618997
##  [705,] 1.0655676
##  [706,] 0.9894063
##  [707,] 0.9681868
##  [708,] 1.0655676
##  [709,] 1.0427419
##  [710,] 0.9681868
##  [711,] 0.9944049
##  [712,] 1.0523772
##  [713,] 1.0085864
##  [714,] 1.0000467
##  [715,] 0.9319773
##  [716,] 1.0575953
##  [717,] 0.9808044
##  [718,] 0.9681868
##  [719,] 1.0655676
##  [720,] 1.0844318
##  [721,] 1.0645854
##  [722,] 1.1618997
##  [723,] 0.9858118
##  [724,] 0.9542853
##  [725,] 0.9894063
##  [726,] 0.9237132
##  [727,] 1.0645854
##  [728,] 0.9674016
##  [729,] 0.9746116
##  [730,] 0.9808044
##  [731,] 0.9894063
##  [732,] 1.0028665
##  [733,] 1.0575953
##  [734,] 0.9674016
##  [735,] 0.9832302
##  [736,] 0.9681868
##  [737,] 1.0575953
##  [738,] 1.1618997
##  [739,] 0.9808044
##  [740,] 0.9976266
##  [741,] 0.9144307
##  [742,] 1.0427419
##  [743,] 0.9746116
##  [744,] 0.9808044
##  [745,] 0.9361535
##  [746,] 1.0758820
##  [747,] 0.9185880
##  [748,] 0.9049066
##  [749,] 0.9681868
##  [750,] 1.0645854
##  [751,] 0.9144307
##  [752,] 0.9858118
##  [753,] 1.1618997
##  [754,] 0.9965271
##  [755,] 0.9319773
##  [756,] 1.1056910
##  [757,] 1.0655676
##  [758,] 1.0645854
##  [759,] 1.0844318
##  [760,] 1.0575953
##  [761,] 0.9185880
##  [762,] 0.9746116
##  [763,] 1.1618997
##  [764,] 1.0028665
##  [765,] 0.9542853
##  [766,] 1.0645854
##  [767,] 1.1006436
##  [768,] 0.9746116
##  [769,] 0.9542853
##  [770,] 1.0844318
##  [771,] 0.9049066
##  [772,] 0.9542853
##  [773,] 0.9894063
##  [774,] 1.0039037
##  [775,] 0.9237132
##  [776,] 1.0575953
##  [777,] 0.9858118
##  [778,] 0.9858118
##  [779,] 1.0655676
##  [780,] 0.9965271
##  [781,] 0.9808044
##  [782,] 1.0427419
##  [783,] 1.0039037
##  [784,] 0.9965271
##  [785,] 0.9319773
##  [786,] 1.0645854
##  [787,] 0.9965271
##  [788,] 1.1006436
##  [789,] 0.9681868
##  [790,] 0.9542853
##  [791,] 1.0844318
##  [792,] 1.0655676
##  [793,] 0.9462140
##  [794,] 1.0655676
##  [795,] 1.0085864
##  [796,] 1.0645854
##  [797,] 0.9361535
##  [798,] 0.9319773
##  [799,] 0.9237132
##  [800,] 1.0028665
##  [801,] 1.0645854
##  [802,] 1.0645854
##  [803,] 1.1006436
##  [804,] 0.9674016
##  [805,] 0.9681868
##  [806,] 0.9681868
##  [807,] 1.0039037
##  [808,] 1.0523772
##  [809,] 0.9746116
##  [810,] 1.0844318
##  [811,] 1.0575953
##  [812,] 1.0000467
##  [813,] 1.1618997
##  [814,] 0.9746116
##  [815,] 1.0427419
##  [816,] 1.0273393
##  [817,] 1.0655676
##  [818,] 0.9361535
##  [819,] 0.9746116
##  [820,] 0.9858118
##  [821,] 0.9144307
##  [822,] 0.9858118
##  [823,] 0.9858118
##  [824,] 1.0427419
##  [825,] 0.9965271
##  [826,] 0.9944049
##  [827,] 0.9144307
##  [828,] 0.9462140
##  [829,] 0.9944049
##  [830,] 0.9361535
##  [831,] 1.0273393
##  [832,] 1.1056910
##  [833,] 0.9681868
##  [834,] 1.0427419
##  [835,] 1.0844318
##  [836,] 0.9237132
##  [837,] 0.9858118
##  [838,] 1.0273393
##  [839,] 0.9832302
##  [840,] 1.0000467
##  [841,] 0.9134597
##  [842,] 1.0427419
##  [843,] 1.0844318
##  [844,] 0.9134597
##  [845,] 1.0758820
##  [846,] 1.0273393
##  [847,] 0.9542853
##  [848,] 0.9542853
##  [849,] 0.9144307
##  [850,] 1.0273393
##  [851,] 0.9542853
##  [852,] 1.0655676
##  [853,] 0.9237132
##  [854,] 0.9049066
##  [855,] 0.9976266
##  [856,] 1.0000467
##  [857,] 1.0028665
##  [858,] 1.0645854
##  [859,] 1.0844318
##  [860,] 1.0844318
##  [861,] 1.1056910
##  [862,] 0.9319773
##  [863,] 0.9944049
##  [864,] 1.0655676
##  [865,] 1.0000467
##  [866,] 0.9237132
##  [867,] 1.0575953
##  [868,] 1.0273393
##  [869,] 1.0028665
##  [870,] 1.1006436
##  [871,] 0.9674016
##  [872,] 1.0427419
##  [873,] 0.9828184
##  [874,] 1.1056910
##  [875,] 1.0000467
##  [876,] 1.0028665
##  [877,] 0.9976266
##  [878,] 1.1056910
##  [879,] 0.9858118
##  [880,] 1.1006436
##  [881,] 1.0523772
##  [882,] 0.9144307
##  [883,] 1.0523772
##  [884,] 1.0427419
##  [885,] 0.9674016
##  [886,] 1.0844318
##  [887,] 0.9858118
##  [888,] 0.9237132
##  [889,] 0.9832302
##  [890,] 1.0085864
##  [891,] 0.9965271
##  [892,] 1.0655676
##  [893,] 1.0273393
##  [894,] 0.9965271
##  [895,] 1.0273393
##  [896,] 0.9828184
##  [897,] 0.9674016
##  [898,] 1.0039037
##  [899,] 0.9542853
##  [900,] 0.9976266
##  [901,] 0.9185880
##  [902,] 1.0039037
##  [903,] 1.0427419
##  [904,] 0.9134597
##  [905,] 0.9319773
##  [906,] 0.9674016
##  [907,] 0.9808044
##  [908,] 1.0028665
##  [909,] 1.0273393
##  [910,] 1.0645854
##  [911,] 1.1618997
##  [912,] 1.0575953
##  [913,] 1.0000467
##  [914,] 0.9832302
##  [915,] 0.9808044
##  [916,] 1.0427419
##  [917,] 0.9858118
##  [918,] 1.0273393
##  [919,] 0.9185880
##  [920,] 0.9976266
##  [921,] 1.0844318
##  [922,] 0.9462140
##  [923,] 0.9134597
##  [924,] 0.9185880
##  [925,] 1.1006436
##  [926,] 1.1618997
##  [927,] 0.9746116
##  [928,] 0.9361535
##  [929,] 0.9049066
##  [930,] 1.0039037
##  [931,] 0.9319773
##  [932,] 0.9049066
##  [933,] 1.0758820
##  [934,] 0.9828184
##  [935,] 0.9185880
##  [936,] 0.9832302
##  [937,] 1.0427419
##  [938,] 0.9237132
##  [939,] 0.9808044
##  [940,] 1.0758820
##  [941,] 1.0575953
##  [942,] 0.9858118
##  [943,] 1.0039037
##  [944,] 0.9134597
##  [945,] 0.9965271
##  [946,] 0.9808044
##  [947,] 0.9049066
##  [948,] 0.9674016
##  [949,] 0.9462140
##  [950,] 0.9319773
##  [951,] 0.9965271
##  [952,] 0.9462140
##  [953,] 1.0655676
##  [954,] 0.9237132
##  [955,] 0.9944049
##  [956,] 1.0000467
##  [957,] 0.9674016
##  [958,] 0.9134597
##  [959,] 0.9944049
##  [960,] 0.9858118
##  [961,] 0.9542853
##  [962,] 0.9462140
##  [963,] 0.9894063
##  [964,] 0.9976266
##  [965,] 0.9361535
##  [966,] 0.9858118
##  [967,] 1.0844318
##  [968,] 0.9185880
##  [969,] 1.0523772
##  [970,] 1.0039037
##  [971,] 1.0645854
##  [972,] 1.0844318
##  [973,] 0.9144307
##  [974,] 1.0575953
##  [975,] 0.9134597
##  [976,] 0.9134597
##  [977,] 0.9808044
##  [978,] 0.9361535
##  [979,] 0.9462140
##  [980,] 1.0273393
##  [981,] 1.1006436
##  [982,] 0.9746116
##  [983,] 0.9976266
##  [984,] 0.9681868
##  [985,] 1.0645854
##  [986,] 0.9049066
##  [987,] 1.0000467
##  [988,] 1.0085864
##  [989,] 0.9858118
##  [990,] 1.0655676
##  [991,] 0.9858118
##  [992,] 0.9134597
##  [993,] 0.9746116
##  [994,] 0.9134597
##  [995,] 0.9681868
##  [996,] 0.9542853
##  [997,] 0.9049066
##  [998,] 0.9542853
##  [999,] 0.9144307
## 
## $model.matrix
##   (Intercept) site1
## 1           1     1
## 2           1    -1
## 3           1     1
## 4           1    -1
## 5           1     1
## 6           1    -1
## 7           1     1
## 8           1    -1
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
```

```r
anova(betadisper(dist_tab_assay,samplesNatSim$site))
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df Sum Sq Mean Sq F value Pr(>F)
## Groups     1   9234  9234.5  1.1703 0.3209
## Residuals  6  47344  7890.6
```

```r
# Exporting results
resOrdered_sp_amb_VS_gm_low_natSim <- sp_amb_VS_gm_low_natSim[order(sp_amb_VS_gm_low_natSim$padj),]
resOrderedDF_sp_amb_VS_gm_low_natSim <- as.data.frame(resOrdered_sp_amb_VS_gm_low_natSim)

write.csv(resOrderedDF_sp_VS_gm, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/juvenile/DESeq2_results_juvenile_sp_VS_gm.csv',sep='\t')
```

```
## Warning in write.csv(resOrderedDF_sp_VS_gm, file = "/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/
## juvenile/DESeq2_results_juvenile_sp_VS_gm.csv", : attempt to set 'sep' ignored
```

```r
write.csv(resOrderedDF_amb_VS_ext, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/juvenile/DESeq2_results_juvenile_amb_VS_ext.csv',sep='\t')
```

```
## Warning in write.csv(resOrderedDF_amb_VS_ext, file = "/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/
## juvenile/DESeq2_results_juvenile_amb_VS_ext.csv", : attempt to set 'sep' ignored
```

```r
write.csv(resOrderedDF_amb_VS_low, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/juvenile/DESeq2_results_juvenile_amb_VS_low.csv',sep='\t')
```

```
## Warning in write.csv(resOrderedDF_amb_VS_low, file = "/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/
## juvenile/DESeq2_results_juvenile_amb_VS_low.csv", : attempt to set 'sep' ignored
```

```r
write.csv(resOrderedDF_low_VS_ext, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/juvenile/DESeq2_results_juvenile_low_VS_ext.csv',sep='\t')
```

```
## Warning in write.csv(resOrderedDF_low_VS_ext, file = "/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/
## juvenile/DESeq2_results_juvenile_low_VS_ext.csv", : attempt to set 'sep' ignored
```

```r
write.csv(resOrderedDF_sp_amb_VS_gm_low_natSim, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/juvenile/DESeq2_results_juvenile_sp_amb_VS_gm_low_natSim.csv',sep='\t')
```

```
## Warning in write.csv(resOrderedDF_sp_amb_VS_gm_low_natSim, file = "/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/
## annotatedGenome/juvenile/DESeq2_results_juvenile_sp_amb_VS_gm_low_natSim.csv", : attempt to set 'sep' ignored
```

```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS 13.2.1
## 
## Matrix products: default
## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] limma_3.50.3                EnhancedVolcano_1.12.0      ashr_2.2-54                 apeglm_1.16.0               tximport_1.22.0            
##  [6] ggvenn_0.1.9                dplyr_1.1.0                 vegan_2.6-2                 lattice_0.20-45             permute_0.9-7              
## [11] gplots_3.1.3                genefilter_1.76.0           RColorBrewer_1.1-3          pheatmap_1.0.12             markdown_1.1               
## [16] ggrepel_0.9.1               ggplot2_3.4.1               BiocManager_1.30.18         devtools_2.4.3              usethis_2.1.6              
## [21] DESeq2_1.34.0               SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0        matrixStats_0.62.0         
## [26] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0              S4Vectors_0.32.4            BiocGenerics_0.40.0        
## [31] knitr_1.39                 
## 
## loaded via a namespace (and not attached):
##   [1] plyr_1.8.7             splines_4.1.2          BiocParallel_1.28.3    digest_0.6.29          invgamma_1.1           htmltools_0.5.2        SQUAREM_2021.1        
##   [8] fansi_1.0.4            magrittr_2.0.3         memoise_2.0.1          cluster_2.1.3          remotes_2.4.2          Biostrings_2.62.0      annotate_1.72.0       
##  [15] extrafont_0.18         extrafontdb_1.0        bdsmatrix_1.3-6        prettyunits_1.1.1      colorspace_2.1-0       blob_1.2.3             xfun_0.31             
##  [22] callr_3.7.0            crayon_1.5.1           RCurl_1.98-1.6         survival_3.3-1         glue_1.6.2             gtable_0.3.1           zlibbioc_1.40.0       
##  [29] XVector_0.34.0         DelayedArray_0.20.0    proj4_1.0-11           pkgbuild_1.3.1         Rttf2pt1_1.3.10        maps_3.4.0             scales_1.2.1          
##  [36] mvtnorm_1.1-3          DBI_1.1.2              Rcpp_1.0.8.3           xtable_1.8-4           emdbook_1.3.12         bit_4.0.4              truncnorm_1.0-8       
##  [43] httr_1.4.3             ellipsis_0.3.2         farver_2.1.1           pkgconfig_2.0.3        XML_3.99-0.9           locfit_1.5-9.5         utf8_1.2.3            
##  [50] tidyselect_1.2.0       labeling_0.4.2         rlang_1.0.6            AnnotationDbi_1.56.2   munsell_0.5.0          tools_4.1.2            cachem_1.0.6          
##  [57] cli_3.6.0              generics_0.1.3         RSQLite_2.2.14         evaluate_0.15          stringr_1.4.0          fastmap_1.1.0          yaml_2.3.5            
##  [64] processx_3.5.3         bit64_4.0.5            fs_1.5.2               caTools_1.18.2         purrr_0.3.4            KEGGREST_1.34.0        nlme_3.1-157          
##  [71] mime_0.12              ash_1.0-15             ggrastr_1.0.1          brio_1.1.3             compiler_4.1.2         rstudioapi_0.13        beeswarm_0.4.0        
##  [78] curl_4.3.2             png_0.1-7              testthat_3.1.4         tibble_3.2.0           geneplotter_1.72.0     stringi_1.7.6          highr_0.9             
##  [85] ps_1.7.0               desc_1.4.1             ggalt_0.4.0            Matrix_1.5-1           vctrs_0.5.2            pillar_1.8.1           lifecycle_1.0.3       
##  [92] bitops_1.0-7           irlba_2.3.5            R6_2.5.1               KernSmooth_2.23-20     vipor_0.4.5            sessioninfo_1.2.2      MASS_7.3-57           
##  [99] gtools_3.9.2.1         pkgload_1.2.4          rprojroot_2.0.3        withr_2.5.0            GenomeInfoDbData_1.2.7 mgcv_1.8-40            parallel_4.1.2        
## [106] coda_0.19-4            rmarkdown_2.14         mixsqp_0.3-43          bbmle_1.0.25           numDeriv_2016.8-1.1    ggbeeswarm_0.6.0
```
