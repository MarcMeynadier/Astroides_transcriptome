# Differential expression on Kallisto data 

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

tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}

# Packages dependance

packageCheckClassic(c('DESeq2','tidyverse','devtools','BiocManager','rhdf5','ggplot2','ggrepel'))
#remotes::install_github("pachterlab/sleuth#260")
#BiocManager::install('tximport', force = TRUE)
BiocManager::install('apeglm')
library('tximport')
library('apeglm')

# Data importation - tximport

samples<-read.table('sample_list.txt',header=T)

files<-paste0(samples$sample,'.tsv')

names(files)<-samples$sample

txi<-tximport(files = files,type='kallisto',txOut=T)

names(txi)

head(txi$counts)

dds<-DESeqDataSetFromTximport(txi,colData=samples,design=~condition)

# Differential expression analysis

dds<-DESeq(dds)
cbind(resultsNames(dds))
res1<-results(dds, name = "Intercept", alpha = 0.05)
res2<-results(dds, name = "condition_pv_pv_pre_vs_gm_gm_pre", alpha = 0.05)
res3<-results(dds, name = "condition_sa_sa_pre_vs_gm_gm_pre", alpha = 0.05) 
summary(res1)
summary(res2)
summary(res3) 

# Exploring the results

#MA-plot
resLFC = lfcShrink(dds, coef = "condition_pv_pv_pre_vs_gm_gm_pre", 
                   type="apeglm")

png("DGE_MA-plot.kallisto.png", width=7, height=5, units = "in", res = 300)
plotMA(resLFC, alpha = 0.05, ylim=c(-6,6), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Principal Component Analysis

# rlog transformation
rld = rlog(dds)

pcaData = plotPCA(rld, intgroup="condition", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-rlog.kallisto.png", width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = condition)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("blue", "red")) +
  geom_text_repel(aes(label = condition), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "rlog transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# vst transformation
vsd = vst(dds)

pcaData = plotPCA(vsd, intgroup=c("individual","paris_classification"), 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png("DGE_PCA-vst.kallisto.png", width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = paris_classification)) + 
  geom_point(size = 2) +theme_bw() + scale_color_manual(values = c("blue", "red")) +
  geom_text_repel(aes(label = individual), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0

p = EnhancedVolcano(data.frame(res), lab = NA, x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "SSA/P vs. Normal",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res), ' variables'),
                    legend=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)

png("DGE_VolcanoPlots.kallisto.png", width=7, height=7, units = "in", res = 300)
print(p)
dev.off()
