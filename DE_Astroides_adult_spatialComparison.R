# Differential expression on Kallisto data 

# Spatial comparison 

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
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
#install_github('cran/heatmap.plus')
library(heatmap.plus)

# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
samplesSingle<-read.table('tximport_design_spatialComparisonSingle.txt',header=T)
samplesPaired<-read.table('tximport_design_spatialComparisonPaired.txt',header=T)
tx2gene<-read.table('tx2gene_fullTranscriptome',header=T)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/5_kallisto/larvaeJuvenileAdultTranscriptome/adult/2_spatialComparison'
outputPath<-paste(scriptPath,'/output/DESeq2/larvaeJuvenileAdultTranscriptome/adult/2_spatialComparison/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
filesPaired<-paste0(samplesPaired$sample,'.tsv')
filesSingle<-paste0(samplesSingle$sample,'.tsv')
names(filesPaired)<-samplesPaired$sample
names(filesSingle)<-samplesSingle$sample
txiPaired<-tximport(files = filesPaired,type='kallisto',tx2gene = tx2gene)
txiSingle<-tximport(files = filesSingle,type='kallisto',tx2gene = tx2gene)
names(txiPaired)
head(txiPaired$counts)
names(txiSingle)
head(txiSingle$counts)
ddsPaired<-DESeqDataSetFromTximport(txiPaired,colData=samplesPaired,design= ~site + experiment)
ddsSingle<-DESeqDataSetFromTximport(txiSingle,colData=samplesSingle,design= ~site + experiment)

# pre-filtering
keepPaired <- rowSums(counts(ddsPaired)) >= 10 
ddsPaired <- ddsPaired[keepPaired,]
keepSingle <- rowSums(counts(ddsSingle)) >= 10 
ddsSingle <- ddsSingle[keepSingle,]

# Differential expression analysis

# Paired end sequences 
ddsPaired<-DESeq(ddsPaired)
cbind(resultsNames(ddsPaired))
res_pv_gm_paired<-results(ddsPaired, contrast=c("site","pv","gm"), alpha = 0.05)
res_sa_gm_paired<-results(ddsPaired, contrast=c("site","sa","gm"), alpha = 0.05)
res_sp_gm_paired<-results(ddsPaired, contrast=c("site","sp","gm"), alpha = 0.05)
res_tro_bck_paired<-results(ddsPaired, contrast=c("experiment","tro","bck"), alpha = 0.05)
summary(res_pv_gm_paired)
summary(res_sa_gm_paired)
summary(res_sp_gm_paired)
summary(res_tro_bck_paired)

# Single end sequences 
ddsSingle<-DESeq(ddsSingle)
cbind(resultsNames(ddsSingle))
res_pv_gm_single<-results(ddsSingle, contrast=c("site","pv","gm"), alpha = 0.05)
res_sp_gm_single<-results(ddsSingle, contrast=c("site","sp","gm"), alpha = 0.05)
res_tro_bck_single<-results(ddsSingle, contrast=c("experiment","tro","bck"), alpha = 0.05)
summary(res_pv_gm_single)
summary(res_sp_gm_single)
summary(res_tro_bck_single)

# Exploring the results

# Results pv VS gm

# MA-plot
# Paired
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonPaired_pv_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsPaired, contrast=c("site","pv","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Single
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonSingle_pv_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsSingle, contrast=c("site","pv","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
# Paired
pCutoff = 0.05
FCcutoff = 1.0
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_pv_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_pv_gm_paired), lab = rownames(data.frame(res_pv_gm_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between pv and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_gm_paired), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Single
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonSingle_pv_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_pv_gm_single), lab = rownames(data.frame(res_pv_gm_single)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between pv and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_gm_single), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results sp VS gm

# MA-plot
# Paired
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonPaired_sp_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsPaired, contrast=c("site","sp","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Single
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonSingle_sp_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsSingle, contrast=c("site","sp","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
# Paired
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_sp_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_sp_gm_paired), lab = rownames(data.frame(res_sp_gm_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_sp_gm_paired), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Single
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonSingle_sp_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_sp_gm_single), lab = rownames(data.frame(res_sp_gm_single)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_sp_gm_single), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results sa VS gm

# MA-plot
# Paired
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonPaired_sa_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsPaired, contrast=c("site","sa","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
# Paired
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_sa_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_sa_gm_paired), lab = rownames(data.frame(res_sa_gm_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sa and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_sa_gm_paired), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results tro vs bck

# MA-plot
# Paired
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonPaired_tro_VS_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsPaired, contrast=c("experiment","tro","bck"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Single
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonSingle_tro_VS_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsSingle, contrast=c("experiment","tro","bck"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes")
dev.off()

# Volcano plot
# Paired
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_tro_VS_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_tro_bck_paired), lab = rownames(data.frame(res_tro_bck_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between gm and sa",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_tro_bck_paired), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Single
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonSingle_tro_VS_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_tro_bck_single), lab = rownames(data.frame(res_tro_bck_single)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between gm and sa",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_tro_bck_single), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()


# Principal Component Analysis

# vst transformation

# Paired
vsdPaired = vst(ddsPaired,blind=T)

pcaData = plotPCA(vsdPaired, intgroup=c("site","experiment"), 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_vst_adult_spatialComparisonPaired.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = site, shape = experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("#ff4040","#6699cc","#9bddff","#000080")) +
  scale_shape_manual(values = c("triangle","circle")) +
  geom_text_repel(aes(label = site), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# Single
vsdSingle = vst(ddsSingle,blind=T)

pcaData = plotPCA(vsdSingle, intgroup=c("site","experiment"), 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_vst_adult_spatialComparisonSingle.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = site, shape = experiment)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("#ff4040","#6699cc","#000080")) +
  scale_shape_manual(values = c("triangle","circle")) +
  geom_text_repel(aes(label = site), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

# heatmap

experiment_colors_paired <- unlist(lapply(samplesPaired$experiment,function(x){
  if(grepl('bck',x)) '#FFC0CB' #pink
  else if(grepl('tro',x)) '#808080' #grey
}))

experiment_colors_single <- unlist(lapply(samplesSingle$experiment,function(x){
  if(grepl('bck',x)) '#FFC0CB' #pink
  else if(grepl('tro',x)) '#808080' #grey
}))

site_colors_paired <- unlist(lapply(samplesPaired$site,function(x){
  if(grepl('gm',x)) '#ff4040' #red
  else if(grepl('pv',x)) '#6699cc' #blue1
  else if(grepl('sa',x)) '#9bddff' #blue2
  else if(grepl('sp',x)) '#000080' #blue3
}))

site_colors_single <- unlist(lapply(samplesSingle$site,function(x){
  if(grepl('gm',x)) '#ff4040' #red
  else if(grepl('pv',x)) '#6699cc' #blue1
  else if(grepl('sp',x)) '#000080' #blue3
}))

myColsPaired <- cbind(experiment_colors_paired,site_colors_paired)
myColsSingle <- cbind(experiment_colors_single,site_colors_single)

# vst transformation

# Paired
topVarGenesVsdPaired <- head(order(rowVars(assay(vsdPaired)), decreasing=TRUE), 50 )
png(paste(outputPath,'DGE_heatmap_adult_spatialComparisonPaired.png',sep=''), width=7, height=7, units = "in", res = 300)
heatmap.3(assay(vsdPaired)[topVarGenesVsdPaired,], trace="none",scale="row",keysize=1,key=T,KeyValueName = "Gene expression",
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7, labCol=F,density.info="none",
          ColSideColors=myColsPaired,xlab="sampling sites & experiment conditions",ylab="genes",margins = c(2,8))
legend(0.95,1,legend=c("gm","pv","sa","sp"),fill=c("#ff4040","#6699cc","#9bddff","#000080"),cex=0.5,xpd=T)
legend(0.95,0.91,legend=c("bck","tro"),fill=c('#FFC0CB','#808080'),cex=0.5,xpd=T)
dev.off()

# Single
topVarGenesVsdSingle <- head(order(rowVars(assay(vsdSingle)), decreasing=TRUE), 50 )
png(paste(outputPath,'DGE_heatmap_vst_adult_spatialComparisonSingle.png',sep=''), width=7, height=7, units = "in", res = 300)
heatmap.3(assay(vsdSingle)[topVarGenesVsdSingle,], trace="none",scale="row",keysize=1,key=T,KeyValueName = "Gene expression",
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7, labCol=F,density.info="none",
          ColSideColors=myColsSingle,xlab="sampling sites & experiment conditions",ylab="genes",margins = c(2,8))
legend(0.95,0.98,legend=c("gm","pv","sp"),fill=c("#ff4040","#6699cc","#000080"),cex=0.5,xpd=T)
legend(0.95,0.91,legend=c("bck","tro"),fill=c('#FFC0CB','#808080'),cex=0.5,xpd=T)
dev.off()

# Exporting results
resOrdered_gm_pv <- res_gm_pv[order(res_gm_pv$pvalue),]
resOrdered_gm_sa <- res_gm_sa[order(res_gm_sa$pvalue),]
head(resOrdered_gm_pv)
head(resOrdered_gm_sa)

resOrderedDF_gm_pv <- as.data.frame(resOrdered_gm_pv)
resOrderedDF_gm_sa <- as.data.frame(resOrdered_gm_sa)
write.csv(resOrderedDF_gm_pv, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_gm_VS_pv.csv',sep=''))
write.csv(resOrderedDF_gm_sa, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_gm_VS_sa.csv',sep=''))

sessionInfo()