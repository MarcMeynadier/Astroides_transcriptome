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

packageCheckClassic(c('DESeq2','devtools','BiocManager','ggplot2','ggrepel','markdown','RColorBrewer','genefilter','gplots','vegan','dplyr'))
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
#BiocManager::install('limma')
#devtools::install_github('cran/GMD')
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library('ggvenn')
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
library('limma')
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
samples<-read.table('tximport_design_juvenile.txt',header=T)
samples2<-read.table('tximport_design_juvenile2.txt',header=T)
samples3<-read.table('tximport_design_juvenile3.txt',header=T)
samplesNatSim<-read.table('tximport_design_juvenile_naturalSimulation.txt',header=T)
tx2gene<-read.table('tx2gene_adultTranscriptome',header=T)
candidateGenes<-read.csv('candidateGenes.csv',header=T,sep=',')
scriptPath <- sub("/[^/]+$", "", scriptPath)
scriptPath <- sub("/[^/]+$", "", scriptPath)
dataPath<-'/data/net/6_kallisto/adultTranscriptome/juvenile'
outputPath<-paste(scriptPath,'/output/DESeq2/adultTranscriptome/juvenile/',sep='')
wdPath<-paste(scriptPath,dataPath,sep='')
setwd(wdPath)

# Data importation - txImport
files<-paste0(samples$samples,'.tsv')
files2<-paste0(samples2$samples,'.tsv')
files3<-paste0(samples3$samples,'.tsv')
filesNatSim<-paste0(samplesNatSim$samples,'.tsv')
names(files)<-samples$samples
names(files2)<-samples2$samples
names(files3)<-samples3$samples
names(filesNatSim)<-samplesNatSim$samples
txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
txi2<-tximport(files = files2,type='kallisto',tx2gene = tx2gene)
txi3<-tximport(files = files3,type='kallisto',tx2gene = tx2gene)
txiNatSim<-tximport(files = filesNatSim,type='kallisto',tx2gene = tx2gene)
names(txi)
names(txi2)
names(txi3)
names(txiNatSim)
head(txi$counts)
head(txi2$counts)
head(txi3$counts)
head(txiNatSim$counts)
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site + pH)
dds2<-DESeqDataSetFromTximport(txi2,colData=samples2,design= ~site + pH)
dds3<-DESeqDataSetFromTximport(txi3,colData=samples3,design= ~site_pH)
ddsNatSim<-DESeqDataSetFromTximport(txiNatSim,colData=samplesNatSim,design= ~site_pH)

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
dds2<-DESeq(dds2)
dds3<-DESeq(dds3)
ddsNatSim<-DESeq(ddsNatSim)
cbind(resultsNames(dds))
cbind(resultsNames(dds2))
cbind(resultsNames(dds3))
cbind(resultsNames(ddsNatSim))
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
sp_amb_VS_gm_low_natSim<-results(ddsNatSim,contrast=c("site_pH","sp_amb","gm_low"),alpha = 0.05)
summary(sp_VS_gm_global)
summary(sp_VS_gm)
summary(amb_VS_ext)
summary(amb_VS_low)
summary(low_VS_ext)
summary(sp_amb_VS_gm_low_natSim)

# Exploring the results

# Results sp VS gm

#MA-plot
png(paste(outputPath,'DGE_MA-plot_juvenile_sp_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(sp_VS_gm,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nsp_VS_gm")
dev.off()
# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
png(paste(outputPath,'DGE_volcanoPlot_juvenile_sp_VS_gm.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(sp_VS_gm), lab = rownames(data.frame(sp_VS_gm)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_VS_gm), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results ext VS amb

#MA-plot
png(paste(outputPath,'DGE_MA-plot_juvenile_amb_VS_ext.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(amb_VS_ext,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\namb_VS_ext")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_juvenile_amb_VS_ext.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(amb_VS_ext), lab = rownames(data.frame(amb_VS_ext)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between ext and amb",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(amb_VS_ext), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()

# Results low VS amb

#MA-plot
png(paste(outputPath,'DGE_MA-plot_juvenile_amb_VS_low.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(amb_VS_low,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\namb_VS_low")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_juvenile_amb_VS_low.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(amb_VS_low), lab = rownames(data.frame(amb_VS_low)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between low and amb",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(amb_VS_low), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()


# Results low VS ext

#MA-plot
png(paste(outputPath,'DGE_MA-plot_juvenile_low_VS_ext.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(low_VS_ext,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nlow_VS_ext")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_juvenile_low_VS_ext.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(low_VS_ext), lab = rownames(data.frame(low_VS_ext)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between low and ext",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(low_VS_ext), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()


# Results natural simulation

#MA-plot
png(paste(outputPath,'DGE_MA-plot_juvenile_natural_simulation.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(sp_amb_VS_gm_low_natSim,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nsp_amb_VS_gm_low")
dev.off()
# Volcano plot
png(paste(outputPath,'DGE_volcanoPlot_juvenile_natural_simulation.png',sep=''), width=7, height=7, units = "in", res = 300)
EnhancedVolcano(data.frame(sp_amb_VS_gm_low_natSim), lab = rownames(data.frame(sp_amb_VS_gm_low_natSim)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp_amb and gm_low",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_amb_VS_gm_low_natSim), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
dev.off()


# Principal Component Analysis
vsd = vst(dds3,blind=T)

pcaData = plotPCA(vsd, intgroup="site_pH", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

pcaData$site_pH = factor(pcaData$site_pH, levels=c("GM_extreme_low","GM_low","GM_ambient","SP_low","SP_ambient"))

png(paste(outputPath,'DGE_PCA_juvenile.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, fill = site_pH)) + 
  geom_point(color="black",pch=21, size=5) + theme_bw() +
  scale_fill_manual(values = c("#D55E00","#E69F00","#FFFF00","#0072B2","#56B4E9")) +
  #ggtitle("Principal Component Analysis of adult corals", subtitle = "may2018 dataset") +
  theme(text = element_text(size=14), legend.position = 'bottom') +
  theme(legend.title=element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()


vsdNatSim = vst(ddsNatSim,blind=T)

pcaData = plotPCA(vsdNatSim, intgroup="site_pH", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

png(paste(outputPath,'DGE_PCA_juvenile_natSim.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = site_pH)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of juvenile corals", subtitle = "Juvenile dataset - Natural conditions simulation") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

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

png(paste(outputPath,'vennDiagramm_juveniles_global.png',sep=''), width=7, height=5, units = "in", res = 300)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#009E73"),
  stroke_size = 0.4, set_name_size = 4
)
dev.off()

resOrdered_sp_VS_gm <- sp_VS_gm[order(sp_VS_gm$padj),]
resOrderedDF_sp_VS_gm <- as.data.frame(resOrdered_sp_VS_gm)
resOrderedDF_sp_VS_gm_venn <- filter(resOrderedDF_sp_VS_gm,padj < 0.05)
resOrderedDF_sp_VS_gm_venn <- list(rownames(resOrderedDF_sp_VS_gm_venn))
resOrderedDF_sp_VS_gm_venn <- unlist(resOrderedDF_sp_VS_gm_venn)

x = list('SP VS GM' = resOrderedDF_sp_VS_gm_venn,'ambient VS low' = resOrderedDF_amb_VS_low_venn)

png(paste(outputPath,'vennDiagramm_juveniles_site_pH.png',sep=''), width=7, height=5, units = "in", res = 300)
ggvenn(
  x, 
  fill_color = c("#8E8E8E", "#E0E0E0"),
  stroke_size = 0.4, set_name_size = 4
)
dev.off()

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

png(paste(outputPath,'vennDiagramm_juveniles_site_amb_VS_low.png',sep=''), width=7, height=5, units = "in", res = 300)
ggvenn(
  x, 
  fill_color = c("#EE4000", "#5CACEE"),
  stroke_size = 0.4, set_name_size = 4
)
dev.off()

x = list('GM ambient VS SP ambient' = resOrderedDF_gm_amb_VS_sp_amb_venn,'GM low VS SP low' = resOrderedDF_gm_low_VS_sp_low_venn)

png(paste(outputPath,'vennDiagramm_juveniles_site_amb_VS_low_2.png',sep=''), width=7, height=5, units = "in", res = 300)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
dev.off()

x = list('GM ambient VS GM extreme low' = resOrderedDF_gm_amb_VS_gm_ext_venn,'GM low VS GM extreme low' = resOrderedDF_gm_low_VS_gm_ext_venn)

png(paste(outputPath,'vennDiagramm_juveniles_gm_ext.png',sep=''), width=7, height=5, units = "in", res = 300)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
dev.off()

# Candidate genes heatmap

# Global
listGenes <- candidateGenes$genes
listGenes2 <- which(rownames(vsd) %in% listGenes)
index <- which(listGenes %in% rownames(vsd))
candidateGenes2 <- candidateGenes[index, ] 
listProt <- candidateGenes2$pfam_annotation
listGenes3 <- candidateGenes2$genes

vsdCandidate <- vsd[listGenes3, ]

labColName <- c('gm_amb','gm_amb','gm_amb','gm_amb','gm_low','gm_low','gm_low','gm_low','gm_ext','gm_ext','gm_ext',
                'gm_ext','sp_amb','sp_amb','sp_amb','sp_amb','sp_low','sp_low','sp_low','sp_low')

colnames(vsdCandidate) <- labColName
rownames(vsdCandidate) <- listProt

topVarGenesVsd <- head(order(rowVars(assay(vsdCandidate)), decreasing=TRUE), 50 )
png(paste(outputPath,'candidateGenes_juveniles_heatmap.png',sep=''), width=7, height=7, units = "in", res = 300)
heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",key.title = "",
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7,density.info="none",
          xlab="sampling sites",ylab="proteins associated to genes",Colv=NA,margins = c(5, 7))

main='Differential expression of 50 most expressed candidates genes\n\nJuveniles'
title(main, cex.main = 0.7)
dev.off()


# Natural conditions simulation
listGenes <- candidateGenes$genes
listGenes2 <- which(rownames(vsdNatSim) %in% listGenes)
index <- which(listGenes %in% rownames(vsdNatSim))
candidateGenes2 <- candidateGenes[index, ] 
listProt <- candidateGenes2$pfam_annotation
listGenes3 <- candidateGenes2$genes

vsdCandidate <- vsdNatSim[listGenes3, ]

labColName <- c('gm_low','gm_low','gm_low','gm_low','sp_amb','sp_amb','sp_amb','sp_amb')

colnames(vsdCandidate) <- labColName
rownames(vsdCandidate) <- listProt

topVarGenesVsd <- head(order(rowVars(assay(vsdCandidate)), decreasing=TRUE), 50 )
png(paste(outputPath,'candidateGenes_juveniles_natSim_heatmap.png',sep=''), width=7, height=7, units = "in", res = 300)
heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",key.title = "",
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7,density.info="none",
          xlab="sampling sites",ylab="proteins associated to genes",Colv=NA,margins = c(5, 7))

main='Differential expression of 50 most expressed candidates genes\n\nJuveniles - Natural conditions simulation'
title(main, cex.main = 0.7)
dev.off()

# Inferences statistics

vsd = vst(dds,blind=T)

count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site + pH, method="euclidian")
anova(betadisper(dist_tab_assay,samples$site))
anova(betadisper(dist_tab_assay,samples$pH))
  
count_tab_assay <- assay(vsdNatSim)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samplesNatSim,dist_tab_assay ~ site_pH, method="euclidian")
anova(betadisper(dist_tab_assay,samplesNatSim$site_pH))

# Exporting results
resOrdered_sp_amb_VS_gm_low_natSim <- sp_amb_VS_gm_low_natSim[order(sp_amb_VS_gm_low_natSim$padj),]
resOrderedDF_sp_amb_VS_gm_low_natSim <- as.data.frame(resOrdered_sp_amb_VS_gm_low_natSim)

write.csv(resOrderedDF_sp_VS_gm, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_sp_VS_gm.csv',sep=''))
write.csv(resOrderedDF_amb_VS_ext, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_amb_VS_ext.csv',sep=''))
write.csv(resOrderedDF_amb_VS_low, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_amb_VS_low.csv',sep=''))
write.csv(resOrderedDF_low_VS_ext, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_low_VS_ext.csv',sep=''))
write.csv(resOrderedDF_sp_amb_VS_gm_low_natSim, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_sp_amb_VS_gm_low_natSim.csv',sep=''))

sessionInfo()
