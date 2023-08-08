# Differential expression on Kallisto data 

# Preliminary samples - 2016 dataset

# Packages and dependence
packageCheckClassic <- function(x){
  # 
  for( i in x ){
    if( ! require( i , character.only = TRUE ) ){
      install.packages( i , dependencies = TRUE )
      require( i , character.only = TRUE )
    }
  }
}

packageCheckClassic(c('DESeq2','adegenet','devtools','BiocManager','ggplot2','ggrepel','markdown','pheatmap','RColorBrewer','genefilter','gplots','vegan','dplyr','limma'))
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("arrayQualityMetrics")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("adegenet")
library('ggvenn')
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
library('BiocManager')
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

# Working environment and data loading
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
candidateGenes<-read.csv('candidateGenes.csv',header=T,sep=',')
samples<-read.table('tximport_design_preliminarySamples.txt',header=T)
dataPath<-'/Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/nov2016'
outputPath<-'/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/output/DESeq2/annotatedGenome/adult/preliminarySamples/'
setwd(dataPath)
data<-list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(data, read.table, skip = 4)
raw_counts <- as.data.frame(sapply(counts.files, function(x) x[ , 4]))
data <- gsub( "Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/nov2016", "", data )
data <- gsub( "_ReadsPerGene.out.tab", "", data )
data <- gsub( "./", "", data )
colnames(raw_counts) <- data
row.names(raw_counts) <- counts.files[[1]]$V1

# DDS object
dds<-DESeqDataSetFromMatrix(countData = raw_counts, colData = samples,design = ~site)

# If data from kallisto
# tx2gene<-read.table('tx2gene_adultTranscriptome',header=T)
# scriptPath <- sub("/[^/]+$", "", scriptPath)
# scriptPath <- sub("/[^/]+$", "", scriptPath)
# dataPath<-'/data/net/6_kallisto/adultTranscriptome/adult/1_preliminarySamples'
# outputPath<-paste(scriptPath,'/output/DESeq2/adultTranscriptome/adult/1_preliminarySamples/',sep='')
# 
# wdPath<-paste(scriptPath,dataPath,sep='')
# setwd(wdPath)
# 
# # Data importation - txImport
# files<-paste0(samples$sample,'.tsv')
# names(files)<-samples$sample
# txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
# names(txi)
# head(txi$counts)
# dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

# Differential expression analysis
dds<-DESeq(dds)
cbind(resultsNames(dds))
res_pv_gm<-results(dds, contrast=c("site","pv","gm"), alpha = 0.05)
res_sa_gm<-results(dds, contrast=c("site","sa","gm"), alpha = 0.05)
res_pv_sa<-results(dds, contrast=c("site","pv","sa"), alpha = 0.05)

# Exploring the results

# Results pv VS gm

#MA-plot
png(paste(outputPath,'DGE_MA-plot_adult_preliminarySamples_pv_VS_gm.png',sep=''), width=7, height=5, units = "in", res = 300)
DESeq2::plotMA(res_pv_gm,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nGarden short : pv VS gm")
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
DESeq2::plotMA(res_sa_gm,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nGarden short : sa VS gm")
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
DESeq2::plotMA(res_pv_sa,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nGarden short : pv VS sa")
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


vsd = vst(dds,blind=T)
mat <- assay(vsd)
mm <- model.matrix(~site,colData(vsd))
mat<-limma::removeBatchEffect(mat,batch1=vsd$site,design=mm)
assay(vsd)<-mat


pcaData = plotPCA(vsd, intgroup="site", 
                  returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))


png(paste(outputPath,'DGE_PCA_adult_preliminarySamples.png',sep=''), width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of adult corals", subtitle = "nov2016 dataset") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
dev.off()

# Venn diagramm 
resOrdered_pv_gm <- res_pv_gm[order(res_pv_gm$padj),]
resOrderedDF_pv_gm <- as.data.frame(resOrdered_pv_gm)
resOrderedDF_pv_gm_venn <- filter(resOrderedDF_pv_gm,padj < 0.05)
resOrderedDF_pv_gm_venn <- list(rownames(resOrderedDF_pv_gm_venn))
resOrderedDF_pv_gm_venn <- unlist(resOrderedDF_pv_gm_venn)

resOrdered_sa_gm <- res_sa_gm[order(res_sa_gm$padj),]
resOrderedDF_sa_gm <- as.data.frame(resOrdered_sa_gm)
resOrderedDF_sa_gm_venn <- filter(resOrderedDF_sa_gm,padj < 0.05)
resOrderedDF_sa_gm_venn <- list(rownames(resOrderedDF_sa_gm_venn))
resOrderedDF_sa_gm_venn <- unlist(resOrderedDF_sa_gm_venn)

resOrdered_pv_sa <- res_pv_sa[order(res_pv_sa$padj),]
resOrderedDF_pv_sa <- as.data.frame(resOrdered_pv_sa)
resOrderedDF_pv_sa_venn <- filter(resOrderedDF_pv_sa,padj < 0.05)
resOrderedDF_pv_sa_venn <- list(rownames(resOrderedDF_pv_sa_venn))
resOrderedDF_pv_sa_venn <- unlist(resOrderedDF_pv_sa_venn)

x = list('pv VS sa' = resOrderedDF_pv_sa_venn, 'sa VS gm' = resOrderedDF_sa_gm_venn, 'pv VS gm' = resOrderedDF_pv_gm_venn)

png(paste(outputPath,'vennDiagramm_preliminarySamples.png',sep=''), width=7, height=5, units = "in", res = 300)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

# Candidate genes heatmap global

listGenes <- candidateGenes$genes
listGenes <- gsub("^([^.]*.[^.]*).*$", "\\1", listGenes)
listGenes <- unique(listGenes)
listGenes2 <- which(rownames(vsd) %in% listGenes)
index <- which(listGenes %in% rownames(vsd))
candidateGenes2 <- candidateGenes[index, ] 
listProt <- candidateGenes2$pfam_annotation
listGenes3 <- candidateGenes2$genes
listGenes3 <- gsub("^([^.]*.[^.]*).*$", "\\1", listGenes3)

removeGenes = c("STRG.35059","STRG.43478")
indexRemoveGenes = which(listGenes3 %in% removeGenes)
listProt <- listProt[-indexRemoveGenes]
listGenes3 <- listGenes3[! listGenes3 %in% c("STRG.35059","STRG.43478")]


vsdCandidate <- vsd[listGenes3, ]

labColName <- samples$originSite_finalSite_experiment
colnames(vsdCandidate) <- labColName
rownames(vsdCandidate) <- paste(listProt,listGenes3,sep=" - ")

topVarGenesVsd <- order(rowVars(assay(vsdCandidate)), decreasing=TRUE)
assayVsdCandidate<-unique(assay(vsdCandidate))
png(paste(outputPath,'candidateGenes_preliminarySamples_heatmap.png',sep=''), width=7, height=7, units = "in", res = 300)
pheatmap(assayVsdCandidate)
dev.off()

png(paste(outputPath,'candidateGenes_preliminarySamples_heatmap_rowScaling.png',sep=''), width=7, height=7, units = "in", res = 300)
pheatmap(assayVsdCandidate,scale = "row")
dev.off()

# Inferences statistics
count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site, method="euclidian")
anova(betadisper(dist_tab_assay,samples$site))

# Exporting results
write.csv(resOrderedDF_pv_gm, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/adult/DESeq2_results_adult_preliminarySamples_pv_VS_gm.csv',sep='\t')
write.csv(resOrderedDF_sa_gm, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/adult/DESeq2_results_adult_preliminarySamples_sa_VS_gm.csv',sep='\t')
write.csv(resOrderedDF_pv_sa, file = '/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/7_deseq2/annotatedGenome/adult/DESeq2_results_adult_preliminarySamples_pv_VS_sa.csv',sep='\t')


heatmap_genes <- function(candidateGenes,vsd,species,genesType,colNames,dataset) {
  
  newColNames = list()
  len = length(newColNames)
  newColNamesSplit = strsplit(sub("^([^_]+_[^_]+)_", "\\1,", colNames), ",")
  c = 1
  for (i in newColNamesSplit) {
    nw = strsplit(sub("^([^_]+_[^_]+_[^_]+)_", "\\1,", i[2]), ",")
    newColNames[[len+c]] = nw[[1]][1]
    c = c + 1
  }
  
  listGenes <- candidateGenes$genes
  listGenes2 <- which(rownames(vsd) %in% listGenes)
  index <- which(listGenes %in% rownames(vsd))
  candidateGenes2 <- candidateGenes[index, ] 
  listProt <- candidateGenes2$pfam_annotation
  listGenes3 <- candidateGenes2$genes
  
  vsdCandidate <- vsd[listGenes3, ]
  
  colnames(vsdCandidate) <- newColNames
  rownames(vsdCandidate) <- listProt
  
  topVarGenesVsd <- head(order(rowVars(assay(vsdCandidate)), decreasing=TRUE), 50 )
  
  heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",
            key.title = "",
            col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7,density.info="none",
            ylab="PFAM annotation of expressed genes",Colv = F,margins = c(6, 7))
  main=paste("Expression level of genes related to",genesType,"in",species,"\n",dataset,sep= " ")
  title(main, cex.main = 0.9)
}


nov2016Cali<-read.csv('/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/8_functionalAnnotation/calicoblasticGenes/nov2016_ortho.csv',header=T,sep=',')
exclusiveGm <- read.csv('/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/data/net/8_functionalAnnotation/calicoblasticGenes/nov2016_gm.csv',header=T,sep=',')

heatmap_genes(nov2016Cali,vsd,'Astroides','calcification',vsd@colData@rownames,"nov2016 orthologs")








sessionInfo()