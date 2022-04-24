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

packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','pheatmap','RColorBrewer','genefilter','gplots','vegan'))
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')

# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
samples<-read.table('tximport_design_preliminarySamples.txt',header=T)
tx2gene<-read.table('tx2gene_adultTranscriptome',header=T)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/5_kallisto/adultTranscriptome/adult/1_preliminarySamples'
outputPath<-paste(scriptPath,'/output/DESeq2/adultTranscriptome/adult/1_preliminarySamples/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
files<-paste0(samples$sample,'.tsv')
names(files)<-samples$sample
txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
names(txi)
head(txi$counts)
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

# Differential expression analysis
dds<-DESeq(dds)
cbind(resultsNames(dds))
res_pv_gm<-results(dds, contrast=c("site","pv","gm"), alpha = 0.05)
res_sa_gm<-results(dds, contrast=c("site","sa","gm"), alpha = 0.05)
res_pv_sa<-results(dds, contrast=c("site","pv","sa"), alpha = 0.05)
summary(res_pv_gm)
summary(res_sa_gm)
summary(res_pv_sa)

# Exploring the results

# Results pv VS gm

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_preliminarySamples_pv_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(dds, contrast=c("site","pv","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), main = "MA-plot for the shrunken log2 fold changes\nPreliminary samples : pv VS gm")
dev.off()

# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
png(paste(outputPath,'DGE_volcanoPlot_adult_preliminarySamples_pv_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_pv_gm), lab = rownames(data.frame(res_pv_gm)), x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "Preliminary samples : pv VS gm",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_gm), ' variables'),
                    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results sa VS gm

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_preliminarySamples_sa_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(dds, contrast=c("site","sa","gm"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes\nPreliminary samples : sa VS gm")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_preliminarySamples_sa_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_sa_gm), lab = rownames(data.frame(res_sa_gm)), x = 'log2FoldChange', y = 'padj',
                    xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                    pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                    title = "Volcano plot", subtitle = "Preliminary samples : sa VS gm",
                    caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_sa_gm), ' variables'),
                    legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                    legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results pv VS sa

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_preliminarySamples_pv_VS_sa.png',sep=''), width=7, height=5, units = "in", res = 300)
resLFC = lfcShrink(dds, contrast=c("site","pv","sa"), 
                   type="ashr")
plotMA(resLFC, alpha = 0.05, ylim=c(-25,25), 
       main = "MA-plot for the shrunken log2 fold changes\nPreliminary samples : pv VS sa")
dev.off()

# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_adult_preliminarySamples_pv_VS_sa.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(res_pv_sa), lab = rownames(data.frame(res_pv_sa)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Preliminary samples : pv VS sa",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(res_pv_sa), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Principal Component Analysis

# vst transformation

vsd = vst(dds,blind=T)

pcaData = plotPCA(vsd, intgroup="site", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_vst_adult_preliminarySamples.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 2) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#6699cc","#9bddff")) +
  geom_text_repel(aes(label = site), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "VST transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  stat_ellipse(level = 0.95)
dev.off()

# heatmap

# vst transformation
topVarGenesVsd <- head(order(rowVars(assay(vsd)), decreasing=TRUE), 50)
png(paste(outputPath,'DGE_heatmap_vst_adult_preliminarySamples.png',sep=''), width=7, height=7, units = "in", res = 300)
heatmap.2(assay(vsd)[topVarGenesVsd,], Rowv=T,trace="none",scale="row",keysize=1, key.par = list(mar=c(3,4,3,0)),
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.5, cexCol=0.7, labCol=F,
          main = "Differentially expressed genes\nin preliminary samples (VST transformation)",
          ColSideColors=c(gm="#ff4040",pv="#6699cc",sa="#9bddff")
          [colData(vsd)$site],xlab="sampling sites", density.info="none",
          ylab="genes",margins=c(2,8))
legend(0.93,1.08,title = "sampling site",legend=c("gm","pv","sa"), 
       fill=c("#ff4040","#6699cc","#9bddff"), cex=0.5, box.lty=1,xpd=T)
dev.off()

# Inferences statistics

count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site, method="euclidian")
anova(betadisper(dist_tab_assay,samples$site))

# Exporting results
resOrdered_pv_gm <- res_pv_gm[order(res_pv_gm$pvalue),]
resOrdered_sa_gm <- res_sa_gm[order(res_sa_gm$pvalue),]
resOrdered_pv_sa <- res_sa_gm[order(res_pv_sa$pvalue),]

head(resOrdered_pv_gm)
head(resOrdered_sa_gm)
head(resOrdered_pv_sa)

resOrderedDF_pv_gm <- as.data.frame(resOrdered_pv_gm)
resOrderedDF_sa_gm <- as.data.frame(resOrdered_sa_gm)
resOrderedDF_pv_sa <- as.data.frame(resOrdered_pv_sa)

write.csv(resOrderedDF_pv_gm, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_pv_VS_gm.csv',sep=''))
write.csv(resOrderedDF_sa_gm, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_sa_VS_gm.csv',sep=''))
write.csv(resOrderedDF_pv_sa, file = paste(scriptPath,'/data/net/6_deseq2/larvaeJuvenileAdultTranscriptome/adult/DESeq2_results_adult_preliminarySamples_pv_VS_sa.csv',sep=''))

sessionInfo()