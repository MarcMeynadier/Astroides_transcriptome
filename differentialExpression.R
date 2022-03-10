# Differential expression on Kallisto data 

# nov2016 dataset - adultTranscriptome

setwd('~/Documents/Projet/code/kallistoResults/adultTranscriptome/adult/nov2016')

# Functions

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

# Packages dependance

packageCheckClassic(c('DESeq2','tidyverse','devtools','BiocManager','rhdf5','ggplot2','ggrepel'))
#remotes::install_github("pachterlab/sleuth#260")
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')

# Data importation - tximport

samples<-read.table('sample_list.txt',header=T)

files<-paste0(samples$sample,'.tsv')

names(files)<-samples$sample

txi<-tximport(files = files,type='kallisto',txOut=T)

names(txi)

head(txi$counts)

dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~ condition)

keep <- rowSums(counts(dds)) >= 10 # pre-filtering
dds <- dds[keep,]

# Differential expression analysis

dds<-DESeq(dds)
cbind(resultsNames(dds))
res_gm_pv<-results(dds, contrast=c("condition","gm","pv"), alpha = 0.05)
res_gm_sa<-results(dds, contrast=c("condition","gm","sa"), alpha = 0.05)
summary(res_gm_pv)
summary(res_gm_sa)

# Exploring the results

# Results gm VS pv

#MA-plot
resLFC = lfcShrink(dds, contrast=c("condition","gm","pv"), 
                   type="ashr")

png("DGE_MA-plot_adult_adultTranscriptome_nov2016_gm_VS_pv.png", width=7, height=5, units = "in", res = 300)
plotMA(resLFC, alpha = 0.05, ylim=c(-6,6), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Principal Component Analysis

# rlog transformation
rld = rlog(dds)

pcaData = plotPCA(rld, intgroup="condition", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-rlog_adult_adultTranscriptome_nov2016_gm_VS_pv.png", width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = condition)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("blue", "red","green")) +
  geom_text_repel(aes(label = condition), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "rlog transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# vst transformation
vsd = vst(dds)

pcaData = plotPCA(vsd, intgroup="condition", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-vst_adult_adultTranscriptome_nov2016_gm_VS_pv.png", width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = condition)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("blue", "red","green")) +
  geom_text_repel(aes(label = condition), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0

p = EnhancedVolcano(data.frame(res_gm_pv), lab = rownames(data.frame(res_gm_pv)), x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "Contrast between gm and pv",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res), ' variables'),
                    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)

png("DGE_VolcanoPlots_adult_adultTranscriptome_nov2016_gm_VS_pv.png", width=7, height=7, units = "in", res = 300)
print(p)
dev.off()

# Results gm VS sa

#MA-plot
resLFC = lfcShrink(dds, contrast=c("condition","gm","sa"), 
                   type="ashr")

png("DGE_MA-plot_adult_adultTranscriptome_nov2016_gm_VS_sa.png", width=7, height=5, units = "in", res = 300)
plotMA(resLFC, alpha = 0.05, ylim=c(-6,6), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Principal Component Analysis

# rlog transformation
rld = rlog(dds)

pcaData = plotPCA(rld, intgroup="condition", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-rlog_adult_adultTranscriptome_nov2016_gm_sa.png", width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = condition)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("blue", "red","green")) +
  geom_text_repel(aes(label = condition), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "rlog transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# vst transformation
vsd = vst(dds)

pcaData = plotPCA(vsd, intgroup="condition", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-vst_adult_adultTranscriptome_nov2016_gm_VS_sa.png", width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = condition)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("blue", "red","green")) +
  geom_text_repel(aes(label = condition), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0

p = EnhancedVolcano(data.frame(res_gm_sa), lab = rownames(data.frame(res)), x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "Contrast between gm and sa",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res), ' variables'),
                    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)

png("DGE_VolcanoPlots_adult_adultTranscriptome_nov2016_gm_VS_sa.png", width=7, height=7, units = "in", res = 300)
print(p)
dev.off()


# Exporting results

resOrdered_gm_pv <- res_gm_pv[order(res$pvalue),]
resOrdered_gm_sa <- res_gm_sa[order(res$pvalue),]
head(resOrdered_gm_pv)
head(resOrdered_gm_sa)

resOrderedDF_gm_pv <- as.data.frame(resOrdered_gm_pv)
resOrderedDF_gm_sa <- as.data.frame(resOrdered_gm_sa)
write.csv(resOrderedDF_gm_pv, file = "DESeq2_results_adult_adultTranscriptome_nov2016_gm_VS_pv.csv")
write.csv(resOrderedDF_gm_sa, file = "DESeq2_results_adult_adultTranscriptome_nov2016_gm_VS_sa.csv")

sessionInfo()
