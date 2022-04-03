# Differential expression on Kallisto data 

# True transplant 

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

packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','RColorBrewer','genefilter','gplots'))
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
#BiocManager::install('limma')
#devtools::install_github('cran/GMD')
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
library('limma')

# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
samples<-read.table('tximport_design_trueTransplant.txt',header=T)
tx2gene<-read.table('tx2gene_fullTranscriptome',header=T)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/5_kallisto/larvaeJuvenileAdultTranscriptome/adult/4_trueTransplant'
outputPath<-paste(scriptPath,'/output/DESeq2/larvaeJuvenileAdultTranscriptome/adult/4_trueTransplant/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
files<-paste0(samples$samples,'.tsv')
names(files)<-samples$samples
txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
names(txi)
head(txi$counts)
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~originSite_finalSite_experiment)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

# Differential expression analysis
dds<-DESeq(dds)
#write.table(counts(dds,normalized=TRUE),file='normalized_count_table_trueTransplant.txt',quote=F)
cbind(resultsNames(dds))
gm_gm_tro_VS_gm_gm_bck<-results(dds, contrast=c("originSite_finalSite_experiment","gm_gm_tro","gm_gm_bck"), alpha = 0.05)
pv_pv_tro_VS_pv_pv_bck<-results(dds, contrast=c("originSite_finalSite_experiment","pv_pv_tro","pv_pv_bck"), alpha = 0.05)
sp_sp_tro_VS_sp_sp_bck<-results(dds, contrast=c("originSite_finalSite_experiment","sp_sp_tro","sp_sp_bck"), alpha = 0.05)
pv_gm_trt_VS_pv_pv_bck<-results(dds, contrast=c("originSite_finalSite_experiment","pv_gm_trt","pv_pv_bck"), alpha = 0.05)
sp_gm_trt_VS_sp_sp_bck<-results(dds, contrast=c("originSite_finalSite_experiment","sp_gm_trt","sp_sp_bck"), alpha = 0.05)
sp_sp_tro_VS_sp_gm_trt<-results(dds, contrast=c("originSite_finalSite_experiment","sp_sp_tro","sp_gm_trt"), alpha = 0.05)
pv_pv_tro_VS_pv_gm_trt<-results(dds, contrast=c("originSite_finalSite_experiment","pv_pv_tro","pv_gm_trt"), alpha = 0.05)
summary(gm_gm_tro_VS_gm_gm_bck)
summary(pv_pv_tro_VS_pv_pv_bck)
summary(sp_sp_tro_VS_sp_sp_bck)
summary(pv_gm_trt_VS_pv_pv_bck)
summary(sp_gm_trt_VS_sp_sp_bck)
summary(sp_sp_tro_VS_sp_gm_trt)
summary(pv_pv_tro_VS_pv_gm_trt)

# Exploring the results

# Results gm_gm_tro VS gm_gm_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_gm_gm_bck_VS_gm_gm_tro.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(gm_gm_tro_VS_gm_gm_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : gm_gm_bck VS gm_gm_tro")
dev.off()
# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_gm_gm_bck_VS_gm_gm_tro.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(gm_gm_tro_VS_gm_gm_bck), lab = rownames(data.frame(gm_gm_tro_VS_gm_gm_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : gm_gm_bck VS gm_gm_tro",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(gm_gm_tro_VS_gm_gm_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results pv_pv_tro VS pv_pv_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_pv_pv_tro_VS_pv_pv_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(pv_pv_tro_VS_pv_pv_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : pv_pv_bck VS pv_pv_tro")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_pv_pv_tro_VS_pv_pv_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(pv_pv_tro_VS_pv_pv_bck), lab = rownames(data.frame(pv_pv_tro_VS_pv_pv_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : pv_pv_bck VS pv_pv_tro",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(pv_pv_tro_VS_pv_pv_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results sp_sp_tro VS sp_sp_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_sp_sp_bck_VS_sp_sp_tro.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(sp_sp_tro_VS_sp_sp_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : sp_sp_bck VS sp_sp_tro")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_sp_sp_bck_VS_sp_sp_tro.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(sp_sp_tro_VS_sp_sp_bck), lab = rownames(data.frame(sp_sp_tro_VS_sp_sp_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : sp_sp_bck VS sp_sp_tro",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_sp_tro_VS_sp_sp_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()


# Results pv_gm_trt VS pv_pv_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_pv_gm_trt_VS_pv_pv_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(pv_gm_trt_VS_pv_pv_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : pv_gm_trt VS pv_pv_bck")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_pv_gm_trt_VS_pv_pv_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(pv_gm_trt_VS_pv_pv_bck), lab = rownames(data.frame(pv_gm_trt_VS_pv_pv_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : pv_gm_trt VS pv_pv_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(pv_gm_trt_VS_pv_pv_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results sp_gm_trt VS sp_sp_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_sp_gm_trt_VS_sp_sp_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(sp_gm_trt_VS_sp_sp_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : sp_gm_trt VS sp_sp_bck")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_sp_gm_trt_VS_sp_sp_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(sp_gm_trt_VS_sp_sp_bck), lab = rownames(data.frame(sp_gm_trt_VS_sp_sp_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : sp_gm_trt VS sp_sp_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_gm_trt_VS_sp_sp_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results pv_pv_tro VS pv_gm_trt

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_pv_pv_tro_VS_pv_gm_trt.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(pv_pv_tro_VS_pv_gm_trt,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : pv_pv_tro VS pv_gm_trt")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_pv_pv_tro_VS_pv_gm_trt.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(pv_pv_tro_VS_pv_gm_trt), lab = rownames(data.frame(pv_pv_tro_VS_pv_gm_trt)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : pv_pv_tro VS pv_gm_trt",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(pv_pv_tro_VS_pv_gm_trt), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results sp_sp_tro VS sp_gm_trt

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_sp_sp_tro_VS_sp_gm_trt.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(sp_sp_tro_VS_sp_gm_trt,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : sp_sp_tro VS sp_gm_trt")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_sp_sp_tro_VS_sp_gm_trt.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(sp_sp_tro_VS_sp_gm_trt), lab = rownames(data.frame(sp_sp_tro_VS_sp_gm_trt)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : sp_sp_tro VS sp_gm_trt",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_sp_tro_VS_sp_gm_trt), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Principal Component Analysis

# vst transformation
vsd = vst(dds,blind=F)

pcaData = plotPCA(vsd, intgroup="originSite_finalSite_experiment", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_adult_trueTransplant.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("#ff0040", "#a40000","#9bddff")) +
  geom_text_repel(aes(label = originSite_finalSite_experiment), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "VST transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# heatmap

# vst transformation
topVarGenesVsd <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 50 )
png(paste(outputPath,'DGE_heatmap_adult_trueTransplant.png',sep=''), width=7, height=7, units = "in", res = 300)
heatmap.2(assay(vsd)[topVarGenesVsd,], Rowv=T,trace="none",scale="row",keysize=1, key.par = list(mar=c(3,4,3,0)),
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.5, cexCol=0.7, labCol=F,
          main = "Differentially expressed genes\nin true transplant experiment",
          ColSideColors=c(gm_pv="#ff0040", gm_sp="#a40000", pv_gm="#6699cc",sp_gm="#9bddff")
          [colData(vsd)$originSite_finalSite],xlab="sampling sites after transplantation", density.info="none",
          ylab="genes",margins=c(2,8))
legend(0.93,1.08,title = "originSite_finalSite",legend=c("gm_pv","gm_sp","pv_gm","sp_gm"), 
       fill=c("#ff0040","#a40000","#6699cc","#9bddff"), cex=0.5, box.lty=1,xpd=T)
dev.off()

# Exporting results
resOrdered_gm_gm_tro_VS_gm_gm_bck <- gm_gm_tro_VS_gm_gm_bck[order(gm_gm_tro_VS_gm_gm_bck$pvalue),]
resOrdered_pv_pv_tro_VS_pv_pv_bck <- pv_pv_tro_VS_pv_pv_bck[order(pv_pv_tro_VS_pv_pv_bck$pvalue),]
resOrdered_sp_sp_tro_VS_sp_sp_bck <- sp_sp_tro_VS_sp_sp_bck[order(sp_sp_tro_VS_sp_sp_bck$pvalue),]
resOrdered_pv_gm_trt_VS_pv_pv_bck <- pv_gm_trt_VS_pv_pv_bck[order(pv_gm_trt_VS_pv_pv_bck$pvalue),]
resOrdered_sp_gm_trt_VS_sp_sp_bck <- sp_gm_trt_VS_sp_sp_bck[order(sp_gm_trt_VS_sp_sp_bck$pvalue),]
resOrdered_sp_sp_tro_VS_sp_gm_trt <- sp_sp_tro_VS_sp_gm_trt[order(sp_sp_tro_VS_sp_gm_trt$pvalue),]
resOrdered_pv_pv_tro_VS_pv_gm_trt <- pv_pv_tro_VS_pv_gm_trt[order(pv_pv_tro_VS_pv_gm_trt$pvalue),]

head(resOrdered_gm_gm_tro_VS_gm_gm_bck)
head(resOrdered_pv_pv_tro_VS_pv_pv_bck)
head(resOrdered_sp_sp_tro_VS_sp_sp_bck)
head(resOrdered_pv_gm_trt_VS_pv_pv_bck)
head(resOrdered_sp_gm_trt_VS_sp_sp_bck)
head(resOrdered_sp_sp_tro_VS_sp_gm_trt)
head(resOrdered_pv_pv_tro_VS_pv_gm_trt)

resOrderedDF_gm_gm_tro_VS_gm_gm_bck <- as.data.frame(resOrdered_gm_gm_tro_VS_gm_gm_bck)
resOrderedDF_pv_pv_tro_VS_pv_pv_bck <- as.data.frame(resOrdered_pv_pv_tro_VS_pv_pv_bck)
resOrderedDF_sp_sp_tro_VS_sp_sp_bck <- as.data.frame(resOrdered_sp_sp_tro_VS_sp_sp_bck)
resOrderedDF_pv_gm_trt_VS_pv_pv_bck <- as.data.frame(resOrdered_pv_gm_trt_VS_pv_pv_bck)
resOrderedDF_sp_gm_trt_VS_sp_sp_bck <- as.data.frame(resOrdered_sp_gm_trt_VS_sp_sp_bck)
resOrderedDF_sp_sp_tro_VS_sp_gm_trt <- as.data.frame(resOrdered_sp_sp_tro_VS_sp_gm_trt)
resOrderedDF_pv_pv_tro_VS_pv_gm_trt <- as.data.frame(resOrdered_pv_pv_tro_VS_pv_gm_trt)

write.csv(resOrderedDF_gm_gm_tro_VS_gm_gm_bck, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_trueTransplant_gm_gm_tro_VS_gm_gm_bck.csv',sep=''))
write.csv(resOrderedDF_pv_pv_tro_VS_pv_pv_bck, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_trueTransplant_pv_pv_tro_VS_pv_pv_bck.csv',sep=''))
write.csv(resOrderedDF_sp_sp_tro_VS_sp_sp_bck, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_trueTransplant_sp_sp_tro_VS_sp_sp_bck.csv',sep=''))
write.csv(resOrderedDF_pv_gm_trt_VS_pv_pv_bck, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_trueTransplant_pv_gm_trt_VS_pv_pv_bck.csv',sep=''))
write.csv(resOrderedDF_sp_gm_trt_VS_sp_sp_bck, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_trueTransplant_sp_gm_trt_VS_sp_sp_bck.csv',sep=''))
write.csv(resOrderedDF_sp_sp_tro_VS_sp_gm_trt, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_trueTransplant_sp_sp_tro_VS_sp_gm_trt.csv',sep=''))
write.csv(resOrderedDF_pv_pv_tro_VS_pv_gm_trt, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_trueTransplant_pv_pv_tro_VS_pv_gm_trt.csv',sep=''))

sessionInfo()
