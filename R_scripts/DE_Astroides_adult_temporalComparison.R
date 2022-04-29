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

packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','RColorBrewer','genefilter','gplots','vegan'))
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
samples <- read.table('tximport_design_temporalComparison.txt',header=T)
samplesSingle<-read.table('tximport_design_temporalComparisonSingle.txt',header=T)
samplesPaired<-read.table('tximport_design_temporalComparisonPaired.txt',header=T)
tx2gene<-read.table('tx2gene_fullTranscriptome',header=T)
scriptPath <- sub("/[^/]+$", "", scriptPath)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/5_kallisto/larvaeJuvenileAdultTranscriptome/adult/3_temporalComparison'
outputPath<-paste(scriptPath,'/output/DESeq2/larvaeJuvenileAdultTranscriptome/adult/3_temporalComparison/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
files<-paste0(samples$sample,'.tsv')
filesPaired<-paste0(samplesPaired$sample,'.tsv')
filesSingle<-paste0(samplesSingle$sample,'.tsv')
names(files)<-samples$samples
names(filesPaired)<-samplesPaired$samples
names(filesSingle)<-samplesSingle$samples
txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
txiPaired<-tximport(files = filesPaired,type='kallisto',tx2gene = tx2gene)
txiSingle<-tximport(files = filesSingle,type='kallisto',tx2gene = tx2gene)
names(txi)
head(txi)
names(txiPaired)
head(txiPaired$counts)
names(txiSingle)
head(txiSingle$counts)
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site + experiment + date)
ddsPaired<-DESeqDataSetFromTximport(txiPaired,colData=samplesPaired,design= ~site + experiment + date)
ddsSingle<-DESeqDataSetFromTximport(txiSingle,colData=samplesSingle,design= ~site + experiment + date)

# pre-filtering
keepPaired <- rowSums(counts(ddsPaired)) >= 10 
ddsPaired <- ddsPaired[keepPaired,]
keepSingle <- rowSums(counts(ddsSingle)) >= 10 
ddsSingle <- ddsSingle[keepSingle,]
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Differential expression analysis

# Paired end sequences 
ddsPaired<-DESeq(ddsPaired)
cbind(resultsNames(ddsPaired))
res_pv_VS_gm_paired<-results(ddsPaired, contrast=c("site","pv","gm"), alpha = 0.05)
res_sa_VS_gm_paired<-results(ddsPaired, contrast=c("site","sa","gm"), alpha = 0.05)
res_sp_VS_gm_paired<-results(ddsPaired, contrast=c("site","sp","gm"), alpha = 0.05)
res_tro_VS_bck_paired<-results(ddsPaired, contrast=c("experiment","tro","bck"), alpha = 0.05)
summary(res_pv_VS_gm_paired)
summary(res_sa_VS_gm_paired)
summary(res_sp_VS_gm_paired)
summary(res_tro_VS_bck_paired)

# Single end sequences 
ddsSingle<-DESeq(ddsSingle)
cbind(resultsNames(ddsSingle))
res_pv_VS_gm_single<-results(ddsSingle, contrast=c("site","pv","gm"), alpha = 0.05)
res_sp_VS_gm_single<-results(ddsSingle, contrast=c("site","sp","gm"), alpha = 0.05)
res_tro_VS_bck_single<-results(ddsSingle, contrast=c("experiment","tro","bck"), alpha = 0.05)
summary(res_pv_VS_gm_single)
summary(res_sp_VS_gm_single)
summary(res_tro_VS_bck_single)

# Both paired end & single end sequences
dds<-DESeq(dds)
cbind(resultsNames(dds))
res_pv_VS_gm<-results(dds, contrast=c("site","pv","gm"), alpha = 0.05)
res_sa_VS_gm<-results(dds, contrast=c("site","sa","gm"), alpha = 0.05)
res_sp_VS_gm<-results(dds, contrast=c("site","sp","gm"), alpha = 0.05)
res_tro_VS_bck<-results(dds, contrast=c("experiment","tro","bck"), alpha = 0.05)
summary(res_pv_VS_gm)
summary(res_sa_VS_gm)
summary(res_sp_VS_gm)
summary(res_tro_VS_bck)

# Exploring the results

# Results pv VS gm

# MA-plot
# Paired
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonPaired_pv_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsPaired, contrast=c("site","pv","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes\nSpatial comparison paired end : pv VS gm")
dev.off()

# Single
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonSingle_pv_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsSingle, contrast=c("site","pv","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes\nSpatial comparison single end : pv VS gm")
dev.off()

# Volcano plot
# Paired
pCutoff = 0.05
FCcutoff = 1.0
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_pv_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_pv_gm_paired), lab = rownames(data.frame(res_pv_gm_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Spatial comparison paired end : pv VS gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_gm_paired), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Single
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonSingle_pv_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_pv_gm_single), lab = rownames(data.frame(res_pv_gm_single)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Spatial comparison single end : pv VS gm",
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
       main = "MA-plot for the shrunken log2 fold changes\nSpatial comparison paired end : sp VS gm")
dev.off()

# Single
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonSingle_sp_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsSingle, contrast=c("site","sp","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes\nSpatial comparison single end : sp VS gm")
dev.off()

# Volcano plot
# Paired
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_sp_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_sp_gm_paired), lab = rownames(data.frame(res_sp_gm_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Spatial comparison paired end : sp VS gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_sp_gm_paired), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Single
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonSingle_sp_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_sp_gm_single), lab = rownames(data.frame(res_sp_gm_single)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Spatial comparison single end : sp VS gm",
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
       main = "MA-plot for the shrunken log2 fold changes\nSpatial comparison paired end : sa VS gm")
dev.off()

# Volcano plot
# Paired
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_sa_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_sa_gm_paired), lab = rownames(data.frame(res_sa_gm_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Spatial comparison paired end : sa VS gm",
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
       main = "MA-plot for the shrunken log2 fold changes\nSpatial comparison paired end : tro VS bck")
dev.off()

# Single
png(paste(outputPath,'DGE_MA-plot_adult_spatialComparisonSingle_tro_VS_bck.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(ddsSingle, contrast=c("experiment","tro","bck"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes\nSpatial comparison single end : tro VS bck")
dev.off()

# Volcano plot
# Paired
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonPaired_tro_VS_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_tro_bck_paired), lab = rownames(data.frame(res_tro_bck_paired)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Spatial comparison paired end : tro VS bck",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_tro_bck_paired), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Single
png(paste(outputPath,'DGE_volcanoPlot_adult_spatialComparisonSingle_tro_VS_bck.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_tro_bck_single), lab = rownames(data.frame(res_tro_bck_single)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Spatial comparison single end : tro VS bck",
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
  geom_text_repel(aes(label = site), nudge_x = -1, nudge_y = 0.2, size = 3, max.overlaps = Inf) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "VST transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  stat_ellipse(level = 0.95)
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
  geom_text_repel(aes(label = site), nudge_x = -1, nudge_y = 0.2, size = 3, max.overlaps = Inf) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "VST transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  stat_ellipse(level = 0.95)
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

# Inferences statistics

vsd = vst(dds,blind=F)
count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site + experiment, method="euclidian")
anova(betadisper(dist_tab_assay,samples$site))
anova(betadisper(dist_tab_assay,samples$experiment))

# Exporting results

resOrdered_pv_VS_gm_paired <- res_pv_VS_gm_paired[order(res_pv_VS_gm_paired$pvalue),]
resOrdered_sa_VS_gm_paired <- res_sa_VS_gm_paired[order(res_sa_VS_gm_paired$pvalue),]
resOrdered_sp_VS_gm_paired <- res_sp_VS_gm_paired[order(res_sp_VS_gm_paired$pvalue),]
resOrdered_tro_VS_bck_paired <- res_tro_VS_bck_paired[order(res_tro_VS_bck_paired$pvalue),]
resOrdered_pv_VS_gm_single <- res_pv_VS_gm_single[order(res_pv_VS_gm_single$pvalue),]
resOrdered_sp_VS_gm_single <- res_sp_VS_gm_single[order(res_sp_VS_gm_single$pvalue),]
resOrdered_tro_VS_bck_single <- res_tro_VS_bck_single[order(res_tro_VS_bck_single$pvalue),]
resOrdered_pv_VS_gm <- res_pv_VS_gm[order(res_pv_VS_gm$pvalue),]
resOrdered_sa_VS_gm <- res_sa_VS_gm[order(res_sa_VS_gm$pvalue),]
resOrdered_sp_VS_gm <- res_sp_VS_gm[order(res_sp_VS_gm$pvalue),]
resOrdered_tro_VS_bck <- res_tro_VS_bck[order(res_tro_VS_bck$pvalue),]

head(resOrdered_pv_VS_gm_paired)
head(resOrdered_sa_VS_gm_paired)
head(resOrdered_sp_VS_gm_paired)
head(resOrdered_tro_VS_bck_paired)
head(resOrdered_pv_VS_gm_single)
head(resOrdered_sp_VS_gm_single)
head(resOrdered_tro_VS_bck_single)
head(resOrdered_pv_VS_gm)
head(resOrdered_sa_VS_gm)
head(resOrdered_sp_VS_gm)
head(resOrdered_tro_VS_bck)

resOrderedDF_pv_VS_gm_paired <- as.data.frame(resOrdered_pv_VS_gm_paired)
resOrderedDF_sa_VS_gm_paired <- as.data.frame(resOrdered_sa_VS_gm_paired)
resOrderedDF_sp_VS_gm_paired <- as.data.frame(resOrdered_sp_VS_gm_paired)
resOrderedDF_tro_VS_bck_paired <- as.data.frame(resOrdered_tro_VS_bck_paired)
resOrderedDF_pv_VS_gm_single <- as.data.frame(resOrdered_pv_VS_gm_single)
resOrderedDF_sp_VS_gm_single <- as.data.frame(resOrdered_sp_VS_gm_single)
resOrderedDF_tro_VS_bck_single <- as.data.frame(resOrdered_tro_VS_bck_single)
resOrderedDF_pv_VS_gm <- as.data.frame(resOrdered_pv_VS_gm)
resOrderedDF_sa_VS_gm <- as.data.frame(resOrdered_sa_VS_gm)
resOrderedDF_sp_VS_gm <- as.data.frame(resOrdered_sp_VS_gm)
resOrderedDF_tro_VS_bck <- as.data.frame(resOrdered_tro_VS_bck)

write.csv(resOrderedDF_pv_VS_gm_paired, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_pv_VS_gm_paired.csv',sep=''))
write.csv(resOrderedDF_sa_VS_gm_paired, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_sa_VS_gm_paired.csv',sep=''))
write.csv(resOrderedDF_sp_VS_gm_paired, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_sp_VS_gm_paired.csv',sep=''))
write.csv(resOrderedDF_tro_VS_bck_paired, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_tro_VS_bck_paired.csv',sep=''))
write.csv(resOrderedDF_pv_VS_gm_single, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_pv_VS_gm_single.csv',sep=''))
write.csv(resOrderedDF_sp_VS_gm_single, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_sp_VS_gm_single.csv',sep=''))
write.csv(resOrderedDF_tro_VS_bck_single, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_tro_VS_bck_single.csv',sep=''))
write.csv(resOrderedDF_pv_VS_gm, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_pv_VS_gm.csv',sep=''))
write.csv(resOrderedDF_sa_VS_gm, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_sa_VS_gm.csv',sep=''))
write.csv(resOrderedDF_sp_VS_gm, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_sp_VS_gm.csv',sep=''))
write.csv(resOrderedDF_tro_VS_bck, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_spatialComparison_tro_VS_bck.csv',sep=''))

sessionInfo()