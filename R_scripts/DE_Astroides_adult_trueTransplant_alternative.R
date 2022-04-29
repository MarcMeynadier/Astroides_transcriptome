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

packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','RColorBrewer','genefilter','gplots','vegan'))
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
samplesBck<-read.table('tximport_design_trueTransplant_bck.txt',header=T)
samplesTro<-read.table('tximport_design_trueTransplant_tro.txt',header=T)
samplesTrt<-read.table('tximport_design_trueTransplant_trt.txt',header=T)
samplesBckTro<-read.table('tximport_design_trueTransplant_tro_bck.txt',header=T)
tx2gene<-read.table('tx2gene_adultTranscriptome',header=T)
scriptPath <- sub("/[^/]+$", "", scriptPath)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/6_kallisto/adultTranscriptome/adult/4_trueTransplant'
outputPath<-paste(scriptPath,'/output/DESeq2/adultTranscriptome/adult/4_trueTransplant/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
filesBck<-paste0(samplesBck$samples,'.tsv')
filesTro<-paste0(samplesTro$samples,'.tsv')
filesTrt<-paste0(samplesTrt$samples,'.tsv')
filesBckTro<-paste0(samplesBckTro$samples,'.tsv')
names(filesBck)<-samplesBck$samples
names(filesTro)<-samplesTro$samples
names(filesTrt)<-samplesTrt$samples
names(filesBckTro)<-samplesBckTro$samples
txiBck<-tximport(files = filesBck,type='kallisto',tx2gene = tx2gene)
txiTro<-tximport(files = filesTro,type='kallisto',tx2gene = tx2gene)
txiTrt<-tximport(files = filesTrt,type='kallisto',tx2gene = tx2gene)
txiBckTro<-tximport(files = filesBckTro,type='kallisto',tx2gene = tx2gene)
ddsBck<-DESeqDataSetFromTximport(txiBck,colData=samplesBck,design= ~originSite_finalSite_experiment)
ddsTro<-DESeqDataSetFromTximport(txiTro,colData=samplesTro,design= ~originSite_finalSite_experiment)
ddsTrt<-DESeqDataSetFromTximport(txiTrt,colData=samplesTrt,design= ~originSite_finalSite_experiment)
ddsBckTro<-DESeqDataSetFromTximport(txiBckTro,colData=samplesBckTro,design= ~originSite_finalSite_experiment)

# pre-filtering
keep <- rowSums(counts(ddsBck)) >= 10 
ddsBck <- ddsBck[keep,]
keep <- rowSums(counts(ddsTro)) >= 10 
ddsTro <- ddsTro[keep,]
keep <- rowSums(counts(ddsTrt)) >= 10 
ddsTrt <- ddsTrt[keep,]
keep <- rowSums(counts(ddsBckTro)) >= 10 
ddsBckTro <- ddsBckTro[keep,]

# Differential expression analysis
ddsBck<-DESeq(ddsBck)
ddsTro<-DESeq(ddsTro)
ddsTrt<-DESeq(ddsTrt)
ddsBckTro<-DESeq(ddsBckTro)
cbind(resultsNames(ddsBck))
cbind(resultsNames(ddsTro))
cbind(resultsNames(ddsTrt))
cbind(resultsNames(ddsBckTro))
pv_pv_bck_VS_gm_gm_bck<-results(ddsBck, contrast=c("originSite_finalSite_experiment","pv_pv_bck","gm_gm_bck"), alpha = 0.05)
sp_sp_bck_VS_gm_gm_bck<-results(ddsBck, contrast=c("originSite_finalSite_experiment","sp_sp_bck","gm_gm_bck"), alpha = 0.05)
pv_pv_bck_VS_sp_sp_bck<-results(ddsBck, contrast=c("originSite_finalSite_experiment","pv_pv_bck","sp_sp_bck"), alpha = 0.05)
pv_pv_tro_VS_gm_gm_tro<-results(ddsTro, contrast=c("originSite_finalSite_experiment","pv_pv_tro","gm_gm_tro"), alpha = 0.05)
sp_sp_tro_VS_gm_gm_tro<-results(ddsTro, contrast=c("originSite_finalSite_experiment","sp_sp_tro","gm_gm_tro"), alpha = 0.05)
pv_pv_tro_VS_sp_sp_tro<-results(ddsTro, contrast=c("originSite_finalSite_experiment","pv_pv_tro","sp_sp_tro"), alpha = 0.05)


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

# Results gm_pv_trt VS pv_pv_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_gm_pv_trt_VS_pv_pv_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(gm_pv_trt_VS_pv_pv_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : gm_pv_trt VS pv_pv_bck")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_gm_pv_trt_VS_pv_pv_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(gm_pv_trt_VS_pv_pv_bck), lab = rownames(data.frame(gm_pv_trt_VS_pv_pv_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : gm_pv_trt VS pv_pv_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(gm_pv_trt_VS_pv_pv_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results gm_sp_trt VS sp_sp_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_gm_sp_trt_VS_sp_sp_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(gm_sp_trt_VS_sp_sp_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : gm_sp_trt VS sp_sp_bck")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_gm_sp_trt_VS_sp_sp_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(gm_sp_trt_VS_sp_sp_bck), lab = rownames(data.frame(gm_sp_trt_VS_sp_sp_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : gm_sp_trt VS sp_sp_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(gm_sp_trt_VS_sp_sp_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results gm_pv_trt VS gm_gm_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_gm_pv_trt_VS_gm_gm_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(gm_pv_trt_VS_gm_gm_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : gm_pv_trt VS gm_gm_bck")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_gm_pv_trt_VS_gm_gm_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(gm_pv_trt_VS_gm_gm_bck), lab = rownames(data.frame(gm_pv_trt_VS_gm_gm_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : gm_pv_trt VS gm_gm_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(gm_pv_trt_VS_gm_gm_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results gm_sp_trt VS gm_gm_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_gm_sp_trt_VS_gm_gm_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(gm_sp_trt_VS_gm_gm_bck,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nTrue transplant : gm_sp_trt VS gm_gm_bck")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_trueTransplant_gm_sp_trt_VS_gm_gm_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(gm_sp_trt_VS_gm_gm_bck), lab = rownames(data.frame(gm_sp_trt_VS_gm_gm_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "True transplant : gm_sp_trt VS gm_gm_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(gm_sp_trt_VS_gm_gm_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Principal Component Analysis

# Background data
vsd = vst(ddsBck,blind=T)

pcaData = plotPCA(vsd, intgroup="originSite_finalSite_experiment", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_adult_trueTransplant_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  #scale_color_manual(values = c("#ff0040", "#a40000","#9bddff")) +
  geom_text_repel(aes(label = originSite_finalSite_experiment), nudge_x = -1, nudge_y = 0.2, size = 3,max.overlaps = Inf) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "VST transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# Transplant origin data
vsd = vst(ddsTro,blind=T)

pcaData = plotPCA(vsd, intgroup="originSite_finalSite_experiment", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_adult_trueTransplant_tro.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  #scale_color_manual(values = c("#ff0040", "#a40000","#9bddff")) +
  geom_text_repel(aes(label = originSite_finalSite_experiment), nudge_x = -1, nudge_y = 0.2, size = 3,max.overlaps = Inf) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "VST transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# Transplant true data
vsd = vst(ddsTrt,blind=T)

pcaData = plotPCA(vsd, intgroup="originSite_finalSite_experiment", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_adult_trueTransplant_trt.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  #scale_color_manual(values = c("#ff0040", "#a40000","#9bddff")) +
  geom_text_repel(aes(label = originSite_finalSite_experiment), nudge_x = -1, nudge_y = 0.2, size = 3,max.overlaps = Inf) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "VST transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()


# Transplant origin and background data
vsd = vst(ddsBckTro,blind=T)

pcaData = plotPCA(vsd, intgroup="originSite_finalSite_experiment", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_adult_trueTransplant_BckTro.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  #scale_color_manual(values = c("#ff0040", "#a40000","#9bddff")) +
  geom_text_repel(aes(label = originSite_finalSite_experiment), nudge_x = -1, nudge_y = 0.2, size = 3,max.overlaps = Inf) +
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
          xlab="sampling sites after transplantation", density.info="none",
          ylab="genes",margins=c(2,8))
#legend(0.93,1.08,title = "originSite_finalSite",legend=c("gm_pv","gm_sp","pv_gm","sp_gm"), 
#       fill=c("#ff0040","#a40000","#6699cc","#9bddff"), cex=0.5, box.lty=1,xpd=T)
dev.off()

# Inferences statistics

count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ originSite_finalSite_experiment, method="euclidian")
anova(betadisper(dist_tab_assay,samples$originSite_finalSite_experiment))

# Exporting results
pv_pv_bck_VS_gm_gm_bck<-results(ddsBck, contrast=c("originSite_finalSite_experiment","pv_pv_bck","gm_gm_bck"), alpha = 0.05)
sp_sp_bck_VS_gm_gm_bck<-results(ddsBck, contrast=c("originSite_finalSite_experiment","sp_sp_bck","gm_gm_bck"), alpha = 0.05)
pv_pv_bck_VS_sp_sp_bck<-results(ddsBck, contrast=c("originSite_finalSite_experiment","pv_pv_bck","sp_sp_bck"), alpha = 0.05)
pv_pv_tro_VS_gm_gm_tro<-results(ddsTro, contrast=c("originSite_finalSite_experiment","pv_pv_tro","gm_gm_tro"), alpha = 0.05)
sp_sp_tro_VS_gm_gm_tro<-results(ddsTro, contrast=c("originSite_finalSite_experiment","sp_sp_tro","gm_gm_tro"), alpha = 0.05)
pv_pv_tro_VS_sp_sp_tro<-results(ddsTro, contrast=c("originSite_finalSite_experiment","pv_pv_tro","sp_sp_tro"), alpha = 0.05)

resOrdered_pv_pv_bck_VS_gm_gm_bck <- pv_pv_bck_VS_gm_gm_bck[order(pv_pv_bck_VS_gm_gm_bck$padj),]
resOrdered_sp_sp_bck_VS_gm_gm_bck <- sp_sp_bck_VS_gm_gm_bck[order(sp_sp_bck_VS_gm_gm_bck$padj),]
resOrdered_pv_pv_bck_VS_sp_sp_bck <- pv_pv_bck_VS_sp_sp_bck[order(pv_pv_bck_VS_sp_sp_bck$padj),]
resOrdered_pv_pv_tro_VS_gm_gm_tro <- pv_pv_tro_VS_gm_gm_tro[order(pv_pv_tro_VS_gm_gm_tro$padj),]
resOrdered_sp_sp_tro_VS_gm_gm_tro <- sp_sp_tro_VS_gm_gm_tro[order(sp_sp_tro_VS_gm_gm_tro$padj),]
resOrdered_pv_pv_tro_VS_sp_sp_tro <- pv_pv_tro_VS_sp_sp_tro[order(pv_pv_tro_VS_sp_sp_tro$padj),]

resOrderedDF_pv_pv_bck_VS_gm_gm_bck <- as.data.frame(resOrdered_pv_pv_bck_VS_gm_gm_bck)
resOrderedDF_sp_sp_bck_VS_gm_gm_bck <- as.data.frame(resOrdered_sp_sp_bck_VS_gm_gm_bck)
resOrderedDF_pv_pv_bck_VS_sp_sp_bck <- as.data.frame(resOrdered_pv_pv_bck_VS_sp_sp_bck)
resOrderedDF_pv_pv_tro_VS_gm_gm_tro <- as.data.frame(resOrdered_pv_pv_tro_VS_gm_gm_tro)
resOrderedDF_sp_sp_tro_VS_gm_gm_tro <- as.data.frame(resOrdered_sp_sp_tro_VS_gm_gm_tro)
resOrderedDF_pv_pv_tro_VS_sp_sp_tro <- as.data.frame(resOrdered_pv_pv_tro_VS_sp_sp_tro)

write.csv(resOrderedDF_pv_pv_bck_VS_gm_gm_bck, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_trueTransplant_pv_pv_bck_VS_gm_gm_bck_subset.csv',sep=''))
write.csv(resOrderedDF_sp_sp_bck_VS_gm_gm_bck, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_trueTransplant_sp_sp_bck_VS_gm_gm_bck_subset.csv',sep=''))
write.csv(resOrderedDF_pv_pv_bck_VS_sp_sp_bck, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_trueTransplant_pv_pv_bck_VS_sp_sp_bck_subset.csv',sep=''))
write.csv(resOrderedDF_pv_pv_tro_VS_gm_gm_tro, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_trueTransplant_pv_pv_tro_VS_gm_gm_tro_subset.csv',sep=''))
write.csv(resOrderedDF_sp_sp_tro_VS_gm_gm_tro, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_trueTransplant_sp_sp_tro_VS_gm_gm_tro_subset.csv',sep=''))
write.csv(resOrderedDF_pv_pv_tro_VS_sp_sp_tro, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/adult/DESeq2_results_adult_trueTransplant_pv_pv_tro_VS_sp_sp_tro_subset.csv',sep=''))

sessionInfo()