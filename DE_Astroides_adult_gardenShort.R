# Differential expression on Kallisto data 

# Garden short 

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
samples<-read.table('tximport_design_gardenShort.txt',header=T)
tx2gene<-read.table('tx2gene',header=T)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/5_kallisto/adult/5_gardenShort'
outputPath<-paste(scriptPath,'/output/DESeq2/5_gardenShort/',sep='')
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
dds$originSite_finalSite_experiment <- relevel(dds$originSite_finalSite_experiment, ref = "gm_gm_bck")
ddsGmRef<-DESeq(dds)
dds$originSite_finalSite_experiment <- relevel(dds$originSite_finalSite_experiment, ref = "pv_pv_bck")
ddsPvRef<-DESeq(dds)
dds$originSite_finalSite_experiment <- relevel(dds$originSite_finalSite_experiment, ref = "sp_sp_bck")
ddsSpRef<-DESeq(dds)
cbind(resultsNames(ddsGmRef))
cbind(resultsNames(ddsPvRef))
cbind(resultsNames(ddsSpRef))
gm_gm_tro_VS_gm_gm_bck<-results(ddsGmRef, contrast=c("originSite_finalSite_experiment","gm_gm_gas","gm_gm_bck"), alpha = 0.05)
pv_gm_trt_VS_pv_pv_bck<-results(ddsPvRef, contrast=c("originSite_finalSite_experiment","pv_gm_trt","pv_pv_bck"), alpha = 0.05)
sp_gm_trt_VS_sp_sp_bck<-results(ddsSpRef, contrast=c("originSite_finalSite_experiment","sp_gm_trt","sp_sp_bck"), alpha = 0.05)
summary(gm_gm_tro_VS_gm_gm_bck)
summary(pv_gm_trt_VS_pv_pv_bck)
summary(sp_gm_trt_VS_sp_sp_bck)

# Exploring the results

# Results gm_gm_tro VS gm_gm_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_gardenShort_gm_gm_bck_VS_gm_gm_tro.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsGmRef, contrast=c("originSite_finalSite_experiment","gm_gm_tro","gm_gm_bck"), 
                   type="ashr")
plotMA(resLFC, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
png(paste(outputPath,'DGE_volcanoPlot_adult_gardenShort_gm_gm_bck_VS_gm_gm_tro.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(gm_gm_tro_VS_gm_gm_bck), lab = rownames(data.frame(gm_gm_tro_VS_gm_gm_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between gm_bck and gm_tro",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(gm_gm_tro_VS_gm_gm_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results pv_gm_trt VS pv_pv_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_gardenShort_pv_gm_trt_VS_pv_pv_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsPvRef, contrast=c("originSite_finalSite_experiment","pv_gm_trt","pv_pv_bck"), 
                   type="ashr")
plotMA(resLFC, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_gardenShort_pv_gm_trt_VS_pv_pv_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(pv_gm_trt_VS_pv_pv_bck), lab = rownames(data.frame(pv_gm_trt_VS_pv_pv_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between pv_gm_trt and pv_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(pv_gm_trt_VS_pv_pv_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results sp_gm_trt VS sp_sp_bck

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_gardenShort_sp_gm_trt_VS_sp_sp_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsSpRef, contrast=c("originSite_finalSite_experiment","sp_gm_trt","sp_sp_bck"), 
                   type="ashr")
plotMA(resLFC, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_gardenShort_sp_gm_trt_VS_sp_sp_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(sp_gm_trt_VS_sp_sp_bck), lab = rownames(data.frame(sp_gm_trt_VS_sp_sp_bck)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp_gm_trt and sp_bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_gm_trt_VS_sp_sp_bck), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()


# Principal Component Analysis

# vst transformation
vsd = vst(ddsGmRef,blind=F)

pcaData = plotPCA(vsd, intgroup="originSite_finalSite_experiment", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_vst_adult_gardenShort.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("#ff0040", "#a40000","#9bddff")) +
  geom_text_repel(aes(label = originSite_finalSite_experiment), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# heatmap

# vst transformation
topVarGenesVsd <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 50 )
png(paste(outputPath,'DGE_heatmap_vst_adult_preliminarySamples.png',sep=''), width=7, height=7, units = "in", res = 300)
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
resOrdered_pv_gm_trt_VS_pv_pv_bck <- pv_gm_trt_VS_pv_pv_bck[order(pv_gm_trt_VS_pv_pv_bck$pvalue),]
resOrdered_sp_gm_trt_VS_sp_sp_bck <- sp_gm_trt_VS_sp_sp_bck[order(sp_gm_trt_VS_sp_sp_bck$pvalue),]

head(resOrdered_gm_gm_tro_VS_gm_gm_bck)
head(resOrdered_pv_gm_trt_VS_pv_pv_bck)
head(resOrdered_sp_gm_trt_VS_sp_sp_bck)

resOrderedDF_gm_gm_tro_VS_gm_gm_bck <- as.data.frame(resOrdered_gm_gm_tro_VS_gm_gm_bck)
resOrderedDF_pv_gm_trt_VS_pv_pv_bck <- as.data.frame(resOrdered_pv_gm_trt_VS_pv_pv_bck)
resOrderedDF_sp_gm_trt_VS_sp_sp_bck <- as.data.frame(resOrdered_sp_gm_trt_VS_sp_sp_bck)

write.csv(resOrdered_gm_gm_tro_VS_gm_gm_bck, file = paste(scriptPath,'/data/net/6_deseq2/adult/DESeq2_results_adult_gardenShort_gm_gm_tro_VS_gm_gm_bck.csv',sep=''))
write.csv(resOrderedDF_pv_gm_trt_VS_pv_pv_bck, file = paste(scriptPath,'/data/net/6_deseq2/adult/DESeq2_results_adult_gardenShort_pv_gm_trt_VS_pv_pv_bck.csv',sep=''))
write.csv(resOrderedDF_sp_gm_trt_VS_sp_sp_bck, file = paste(scriptPath,'/data/net/6_deseq2/adult/DESeq2_results_adult_gardenShort_sp_gm_trt_VS_sp_sp_bck.csv',sep=''))

sessionInfo()