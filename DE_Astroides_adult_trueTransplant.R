# Differential expression on Kallisto data 

# Preliminary samples 

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
library('GMD')

# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
samples<-read.table('tximport_design_trueTransplant.txt',header=T)
tx2gene<-read.table('tx2gene',header=T)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/5_kallisto/adult/4_trueTransplant'
outputPath<-paste(scriptPath,'/output/DESeq2/4_trueTransplant/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
files<-paste0(samples$samples,'.tsv')
names(files)<-samples$samples
txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
names(txi)
head(txi$counts)
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~originSite_finalSite)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

# Differential expression analysis
dds<-DESeq(dds)
cbind(resultsNames(dds))
res_gm_sp_VS_gm_pv<-results(dds, name="originSite_finalSite_gm_sp_vs_gm_pv", alpha = 0.05)
res_pv_gm_VS_gm_pv<-results(dds, name="originSite_finalSite_pv_gm_vs_gm_pv", alpha = 0.05)
res_sp_gm_VS_gm_pv<-results(dds, name="originSite_finalSite_sp_gm_vs_gm_pv", alpha = 0.05)
summary(res_gm_sp_VS_gm_pv)
summary(res_pv_gm_VS_gm_pv)
summary(res_sp_gm_VS_gm_pv)

# Exploring the results

# Results gm_sp VS gm_pv

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_trueTransplant_gm_VS_pv.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(dds, coef="originSite_finalSite_pv_gm_vs_gm_pv", 
                   type="apeglm")
plotMA(resLFC, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
png(paste(outputPath,'DGE_volcanoPlot_adult_preliminarySamples_gm_VS_pv.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_pv_gm_VS_gm_pv), lab = rownames(data.frame(res_pv_gm_VS_gm_pv)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between gm and pv",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_gm_VS_gm_pv), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results gm VS sa

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_preliminarySamples_gm_VS_sa.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(dds, contrast=c("site","gm","sa"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_preliminarySamples_gm_VS_sa.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_gm_sa), lab = rownames(data.frame(res_gm_sa)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between gm and sa",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_gm_sa), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Principal Component Analysis

# rlog transformation
rld = rlog(dds)

pcaData = plotPCA(rld, intgroup="originSite_finalSite", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_rlog_adult_preliminarySamples.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("blue", "red","green","yellow")) +
  geom_text_repel(aes(label = originSite_finalSite), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "rlog transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# vst transformation
vsd = vst(dds)

pcaData = plotPCA(vsd, intgroup="originSite_finalSite", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_vst_adult_preliminarySamples.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("blue", "red","green","yellow")) +
  geom_text_repel(aes(label = originSite_finalSite), nudge_x = -1, nudge_y = 0.2, size = 3) +
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
          main = "Differentially expressed genes\nin true transplant experiment (vst transformation)",
          ColSideColors=c(gm_pv="#ff4040", gm_sp="green", pv_gm="yellow",sp_gm="purple")
          [colData(vsd)$originSite_finalSite],xlab="sampling sites", density.info="none",
          ylab="genes",margins=c(2,8))
legend(0.93,1.02,title = "originSite_finalSite",legend=c("gm_pv","gm_sp","pv_gm","sp_gm"), 
       fill=c("#ff4040","green","yellow","purple"), cex=0.5, box.lty=1,xpd=T)
dev.off()

# Exporting results
resOrdered_gm_pv <- res_gm_pv[order(res_gm_pv$pvalue),]
resOrdered_gm_sa <- res_gm_sa[order(res_gm_sa$pvalue),]
head(resOrdered_gm_pv)
head(resOrdered_gm_sa)

resOrderedDF_gm_pv <- as.data.frame(resOrdered_gm_pv)
resOrderedDF_gm_sa <- as.data.frame(resOrdered_gm_sa)
write.csv(resOrderedDF_gm_pv, file = paste(scriptPath,'/data/net/6_deseq2/adult/1_preliminarySamples/DESeq2_results_adult_preliminarySamples_gm_VS_pv.csv',sep=''))
write.csv(resOrderedDF_gm_sa, file = paste(scriptPath,'/data/net/6_deseq2/adult/1_preliminarySamples/DESeq2_results_adult_preliminarySamples_gm_VS_sa.csv',sep=''))

sessionInfo()