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

packageCheckClassic(c('gridExtra','DESeq2','adegenet','devtools','BiocManager','ggplot2','ggrepel','markdown','pheatmap','RColorBrewer','genefilter','gplots','vegan','dplyr','limma'))
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
library('pairwiseAdonis')
library('BiocManager')
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")


# Functions

trimFunction <- function(resLFC,p_adj_cut,lfc_cut) {
  resOrdered<-resLFC[order(resLFC$padj),]
  y<-na.omit(resOrdered)
  resTrim<-y[y$padj < p_adj_cut,]
  resTrim <- resTrim[ which( resTrim$log2FoldChange > lfc_cut | resTrim$log2FoldChange < -lfc_cut), ]
  resTrim <- as.data.frame(resTrim)
  resTrim$genes <- rownames(resTrim)
  return(resTrim)
}

applyAnnot <- function(inputDf,inputDfAnnot) {
  gene_and_annot_col <- character(0)
  inputDf_annot_genes <- inputDfAnnot$genes
  inputDf_annot_pfam <- inputDfAnnot$pfam_annotation
  for (i in 1:length(inputDf_annot_genes)) {
    gene_and_annot_col <- c(gene_and_annot_col, 
                            paste(inputDf_annot_genes[i], " - ", inputDf_annot_pfam[i]))
  }
  inputDf_genes <- inputDf$genes
  new_gene_annot_col <- character(0)
  for (i in inputDf_genes) {
    gene <- i
    for (j in gene_and_annot_col) {
      j_split <- strsplit(j, " ", fixed = TRUE)[[1]]
      if (i %in% j_split[1]) {
        gene <- j
      }
    }
    new_gene_annot_col <- c(new_gene_annot_col, gene)
  }
  inputDf$genes <- new_gene_annot_col
  return(inputDf)
}

heatmapFunction <- function(newColNames,commonGenes,commonGenesAll,vsd,nameCode) {
  vsdCommonGm <- vsd[commonGenes$genes,]
  vsdCommonGm@colData@rownames = newColNames
  annot_genes = commonGenesAll$genes
  genes = names(vsdCommonGm)
  
  new_names = list()
  for (i in annot_genes) {
    gene = i
    for (j in genes) {
      if (i == j) {
        gene = j
      }
    }
    new_names <- c(new_names,gene)
  }
  names(vsdCommonGm) <- new_names
  pheatmap(assay(vsdCommonGm),main=nameCode,scale="row", cluster_rows=T, show_rownames=T,
           cluster_cols=FALSE,cellwidth = 20,fontsize_row = 4)
  return(vsdCommonGm)
}

similarityFunction <- function(dataset1,dataset2) {
  c = 0
  for (i in rownames(dataset1)) {
    if (i %in% rownames(dataset2)) {
      c = c + 1
    }
  }
  sim = (c / as.integer(length(rownames(dataset2)))) * 100
  return(sim)
}


# Nov2016 dataset

# Working environment and data loading
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
A_pfam<-read.csv('Astroides_pfam.csv',header=T,sep=',')
samples<-read.table('tximport_design_preliminarySamples.txt',header=T)
dataPath<-'/Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/nov2016'
outputPath<-'/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/output/DESeq2/annotatedGenome/adult/preliminarySamples/'
setwd(dataPath)
data<-list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(data, read.table, skip = 4)
raw_counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2]))
data <- gsub( "Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/nov2016", "", data )
data <- gsub( "_ReadsPerGene.out.tab", "", data )
data <- gsub( "./", "", data )
colnames(raw_counts) <- data
row.names(raw_counts) <- counts.files[[1]]$V1

# DDS object
dds<-DESeqDataSetFromMatrix(countData = raw_counts, colData = samples,design = ~site)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

# Differential expression analysis
dds<-DESeq(dds)
cbind(resultsNames(dds))


# General 

vsdNov2016 = vst(dds,blind=T)

pcaData = plotPCA(vsdNov2016, intgroup="site",returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 12) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Home environments - November 2016") +
  theme(text = element_text(size=24),legend.text = element_text(size=20), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 



# Common gm

res_pv_gm_lfc2_intra<-results(dds, contrast=c("site","pv","gm"), lfcThreshold = 1,alpha = 0.05)
res_pv_gm_lfc2_intra<-trimFunction(res_pv_gm_lfc2_intra,0.05,1)

res_sa_gm_lfc2_intra<-results(dds, contrast=c("site","sa","gm"), lfcThreshold = 1,alpha = 0.05)
res_sa_gm_lfc2_intra<-trimFunction(res_sa_gm_lfc2_intra,0.05,1)

common_gm_lfc2_intra <- merge(res_pv_gm_lfc2_intra,res_sa_gm_lfc2_intra,by="genes")
common_gm_annot_lfc2_intra = merge(common_gm_lfc2_intra,A_pfam,by="genes")
common_gm_all_lfc2_intra <- applyAnnot(common_gm_lfc2_intra,common_gm_annot_lfc2_intra)

newColNames = c("gm","gm","gm","gm","pv","pv","pv","sp","sp","sp")
vsd_common_gm_lfc2_intra_nov2016 = heatmapFunction(newColNames,common_gm_lfc2_intra,common_gm_all_lfc2_intra,vsdNov2016,"Nov2016 common GM genes - Intra method")


pcaData = plotPCA(vsd_common_gm_lfc2_intra_nov2016, intgroup="site",returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of adult corals", subtitle = "nov2016 dataset - Common GM DEGs") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# Common sa
res_pv_sa_lfc2_intra<-results(dds, contrast=c("site","pv","sa"), lfcThreshold = 1,alpha = 0.05)
res_pv_sa_lfc2_intra<-trimFunction(res_pv_sa_lfc2_intra,0.05,1)

res_gm_sa_lfc2_intra<-results(dds, contrast=c("site","gm","sa"), lfcThreshold = 1,alpha = 0.05)
res_gm_sa_lfc2_intra<-trimFunction(res_gm_sa_lfc2_intra,0.05,1)

common_sa_lfc2_intra <- merge(res_pv_sa_lfc2_intra,res_gm_sa_lfc2_intra,by="genes")
common_sa_annot_lfc2_intra = merge(common_sa_lfc2_intra,A_pfam,by="genes")
common_sa_all_lfc2_intra <- applyAnnot(common_sa_lfc2_intra,common_sa_annot_lfc2_intra)

newColNames = c("gm","gm","gm","gm","pv","pv","pv","sp","sp","sp")
vsd_common_sa_lfc2_intra_nov2016 = heatmapFunction(newColNames,common_sa_lfc2_intra,common_sa_all_lfc2_intra,vsdNov2016,"Nov2016 common SA genes - Intra method")


pcaData = plotPCA(vsd_common_sa_lfc2_intra_nov2016, intgroup="site",returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of adult corals", subtitle = "nov2016 dataset - Common GM DEGs") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# Common pv
res_sa_pv_lfc2_intra<-results(dds, contrast=c("site","sa","pv"), lfcThreshold = 1,alpha = 0.05)
res_sa_pv_lfc2_intra<-trimFunction(res_pv_gm_lfc2_intra,0.05,1)

res_gm_pv_lfc2_intra<-results(dds, contrast=c("site","gm","pv"), lfcThreshold = 1,alpha = 0.05)
res_gm_pv_lfc2_intra<-trimFunction(res_gm_pv_lfc2_intra,0.05,1)

common_pv_lfc2_intra <- merge(res_sa_pv_lfc2_intra,res_gm_pv_lfc2_intra,by="genes")
common_pv_annot_lfc2_intra = merge(common_pv_lfc2_intra,A_pfam,by="genes")
common_pv_all_lfc2_intra <- applyAnnot(common_pv_lfc2_intra,common_pv_annot_lfc2_intra)

newColNames = c("gm","gm","gm","gm","pv","pv","pv","sa","sa","sa")
vsd_common_pv_lfc2_intra_nov2016 = heatmapFunction(newColNames,common_pv_lfc2_intra,common_pv_all_lfc2_intra,vsdNov2016,"Nov2016 common PV genes - Intra method")


pcaData = plotPCA(vsd_common_pv_lfc2_intra_nov2016, intgroup="site",returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of adult corals", subtitle = "nov2016 dataset - Common GM DEGs") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

# Inferences statistics
count_tab_assay <- assay(vsdNov2016)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samples,dist_tab_assay ~ site, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samples$site)
pair.mod
mod <- betadisper(dist_tab_assay,samples$site)
mod.HSD <- TukeyHSD(mod)
mod.HSD


# May 2018 dataset

# Working environment and data loading
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
A_pfam<-read.csv('Astroides_pfam.csv',header=T,sep=',')
samplesBck<-read.table('tximport_design_trueTransplant_bck.txt',header=T)
samplesTroBck<-read.table('tximport_design_trueTransplant_tro_bck.txt',header=T)
samplesTroBck2<-read.table('tximport_design_trueTransplant_tro_bck_2.txt',header=T)
dataPath<-'/Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/may2018'
outputPath<-'/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/output/DESeq2/annotatedGenome/adult/trueTransplant/'
setwd(dataPath)
data<-list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(data, read.table, skip = 4)
raw_counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2]))
data <- gsub( "Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/may2018", "", data )
data <- gsub( "_ReadsPerGene.out.tab", "", data )
data <- gsub( "./", "", data )
colnames(raw_counts) <- data
row.names(raw_counts) <- counts.files[[1]]$V1
raw_counts_bck <- raw_counts[,grep("bck", colnames(raw_counts))]
raw_counts_bck_tro <- raw_counts[,grep("bck|tro", colnames(raw_counts))] 

# DDS object
ddsBck<-DESeqDataSetFromMatrix(countData = raw_counts_bck, colData = samplesBck,design = ~site)
ddsBckTro<-DESeqDataSetFromMatrix(countData = raw_counts_bck_tro, colData = samplesTroBck,design = ~site + experiment)
ddsBckTro2<-DESeqDataSetFromMatrix(countData = raw_counts_bck_tro, colData = samplesTroBck2,design = ~originSite_finalSite_experiment) 

# pre-filtering
keep <- rowSums(counts(ddsBck)) >= 10 
ddsBck <- ddsBck[keep,]
keep <- rowSums(counts(ddsBckTro)) >= 10 
ddsBckTro <- ddsBckTro[keep,]
keep <- rowSums(counts(ddsBckTro2)) >= 10 
ddsBckTro2 <- ddsBckTro2[keep,]

# Differential expression analysis
ddsBck<-DESeq(ddsBck)
cbind(resultsNames(ddsBck))
ddsBckTro<-DESeq(ddsBckTro)
cbind(resultsNames(ddsBckTro))
ddsBckTro2<-DESeq(ddsBckTro2)
cbind(resultsNames(ddsBckTro2))


# General 
vsdMay2018Bck = vst(ddsBck,blind=T)
vsdMay2018BckTro = vst(ddsBckTro,blind=T)
vsdMay2018BckTro2 = vst(ddsBckTro2,blind=T)

pcaData = plotPCA(vsdMay2018Bck, intgroup="site",returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = site)) + 
  geom_point(size = 12) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Home environments - May 2018") +
  theme(text = element_text(size=24),legend.text = element_text(size=20), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

pcaData = plotPCA(vsdMay2018BckTro, intgroup=c("site","experiment"),returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = site,shape = experiment)) + 
  geom_point(size = 12) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  scale_shape_manual(values = c("circle", "triangle")) +
  geom_point() +
  ggtitle("Home environments & transplant origin - May 2018") +
  theme(text = element_text(size=24),legend.text = element_text(size=20), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# Common gm
res_pv_gm_lfc2<-results(ddsBckTro2, contrast=c("originSite_finalSite_experiment","pv_pv_bck","gm_gm_bck"), lfcThreshold = 1,alpha = 0.05)
res_pv_gm_lfc2<-trimFunction(res_pv_gm_lfc2,0.05,1)

res_sp_gm_lfc2<-results(ddsBckTro2, contrast=c("originSite_finalSite_experiment","sp_sp_bck","gm_gm_bck"), lfcThreshold = 1,alpha = 0.05)
res_sp_gm_lfc2<-trimFunction(res_sp_gm_lfc2,0.05,1)

common_gm_lfc2 <- merge(res_pv_gm_lfc2,res_sp_gm_lfc2,by="genes")
common_gm_annot_lfc2 = merge(common_gm_lfc2,A_pfam,by="genes")
common_gm_all_lfc2 <- applyAnnot(common_gm_lfc2,common_gm_annot_lfc2)

newColNames = c("gm_gm_bck","gm_gm_bck","gm_gm_bck","gm_gm_tro","gm_gm_tro","gm_gm_tro","gm_gm_tro","gm_gm_tro","pv_pv_bck","pv_pv_bck","pv_pv_bck","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","sp_sp_bck","sp_sp_bck","sp_sp_bck","sp_sp_tro","sp_sp_tro","sp_sp_tro","sp_sp_tro","sp_sp_tro")
vsd_common_gm_may2018 = heatmapFunction(newColNames,common_gm_lfc2,common_gm_all_lfc2,vsdMay2018BckTro2,"May2018 common GM genes")


# Common pv
res_gm_pv_lfc2<-results(ddsBckTro2, contrast=c("originSite_finalSite_experiment","gm_gm_bck","pv_pv_bck"), lfcThreshold = 1,alpha = 0.05)
res_gm_pv_lfc2<-trimFunction(res_gm_pv_lfc2,0.05,1)

res_sp_pv_lfc2<-results(ddsBckTro2, contrast=c("originSite_finalSite_experiment","sp_sp_bck","pv_pv_bck"), lfcThreshold = 1,alpha = 0.05)
res_sp_pv_lfc2<-trimFunction(res_sp_pv_lfc2,0.05,1)

common_pv_lfc2 <- merge(res_gm_pv_lfc2,res_sp_pv_lfc2,by="genes")
common_pv_annot_lfc2 = merge(common_pv_lfc2,A_pfam,by="genes")
common_pv_all_lfc2 <- applyAnnot(common_pv_lfc2,common_pv_annot_lfc2)

newColNames = c("gm_gm_bck","gm_gm_bck","gm_gm_bck","gm_gm_tro","gm_gm_tro","gm_gm_tro","gm_gm_tro","gm_gm_tro","pv_pv_bck","pv_pv_bck","pv_pv_bck","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","sp_sp_bck","sp_sp_bck","sp_sp_bck","sp_sp_tro","sp_sp_tro","sp_sp_tro","sp_sp_tro","sp_sp_tro")
vsd_common_pv_may2018 = heatmapFunction(newColNames,common_pv_lfc2,common_pv_all_lfc2,vsdMay2018BckTro2,"May2018 common PV genes")


# Common sp
res_gm_sp_lfc2<-results(ddsBckTro2, contrast=c("originSite_finalSite_experiment","gm_gm_bck","sp_sp_bck"), lfcThreshold = 1,alpha = 0.05)
res_gm_sp_lfc2<-trimFunction(res_gm_sp_lfc2,0.05,1)

res_pv_sp_lfc2<-results(ddsBckTro2, contrast=c("originSite_finalSite_experiment","pv_pv_bck","sp_sp_bck"), lfcThreshold = 1,alpha = 0.05)
res_pv_sp_lfc2<-trimFunction(res_pv_sp_lfc2,0.05,1)

common_sp_lfc2 <- merge(res_gm_sp_lfc2,res_pv_sp_lfc2,by="genes")
common_sp_annot_lfc2 = merge(common_sp_lfc2,A_pfam,by="genes")
common_sp_all_lfc2 <- applyAnnot(common_sp_lfc2,common_sp_annot_lfc2)

newColNames = c("gm_gm_bck","gm_gm_bck","gm_gm_bck","gm_gm_tro","gm_gm_tro","gm_gm_tro","gm_gm_tro","gm_gm_tro","pv_pv_bck","pv_pv_bck","pv_pv_bck","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","pv_pv_tro","sp_sp_bck","sp_sp_bck","sp_sp_bck","sp_sp_tro","sp_sp_tro","sp_sp_tro","sp_sp_tro","sp_sp_tro")
vsd_common_sp_may2018 = heatmapFunction(newColNames,common_sp_lfc2,common_sp_all_lfc2,vsdMay2018BckTro2,"May2018 common SP genes")


# Inferences statistics
count_tab_assay <- assay(vsdMay2018Bck)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samplesBck,dist_tab_assay ~ site, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesBck$site)
pair.mod
mod <- betadisper(dist_tab_assay,samplesBck$site)
mod.HSD <- TukeyHSD(mod)
mod.HSD


count_tab_assay <- assay(vsdMay2018BckTro)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samplesTroBck,dist_tab_assay ~ site + experiment, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesTroBck2$originSite_finalSite_experiment)
pair.mod
mod <- betadisper(dist_tab_assay,samplesTroBck2$originSite_finalSite_experiment)
mod.HSD <- TukeyHSD(mod)
mod.HSD


count_tab_assay <- assay(vsdMay2018BckTro2)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samplesTroBck2,dist_tab_assay ~ originSite_finalSite_experiment, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesTroBck2$originSite_finalSite_experiment)
pair.mod
mod <- betadisper(dist_tab_assay,samplesTroBck2$originSite_finalSite_experiment)
mod.HSD <- TukeyHSD(mod)
mod.HSD




# September 2018 dataset

# Working environment and data loading
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
A_pfam<-read.csv('Astroides_pfam.csv',header=T,sep=',')
samplesBck<-read.table('tximport_design_gardenShort_bck.txt',header=T)
samplesGsBck<-read.table('tximport_design_gardenShort_gs_bck.txt',header=T)
samplesGsBck2<-read.table('tximport_design_gardenShort_gs_bck_2.txt',header=T)
dataPath<-'/Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/sept2018'
outputPath<-'/Users/mmeynadier/Documents/Astroides/comparative_transcriptomics_astroides/output/DESeq2/annotatedGenome/adult/gardenShort/'
setwd(dataPath)
data<-list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(data, read.table, skip = 4)
raw_counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2]))
data <- gsub( "Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/sept2018", "", data )
data <- gsub( "_ReadsPerGene.out.tab", "", data )
data <- gsub( "./", "", data )
colnames(raw_counts) <- data
row.names(raw_counts) <- counts.files[[1]]$V1
raw_counts_bck <- raw_counts[,grep("bck", colnames(raw_counts))]
raw_counts_bck_gs <- raw_counts[,grep("bck|gm_gm_gas|pv_pv_gas|sp_sp_gas", colnames(raw_counts))] 

# DDS object
ddsBck<-DESeqDataSetFromMatrix(countData = raw_counts_bck, colData = samplesBck,design = ~originSite_finalSite_experiment)
ddsBckGs<-DESeqDataSetFromMatrix(countData = raw_counts_bck_gs, colData = samplesGsBck,design = ~site + experiment)
ddsBckGs2<-DESeqDataSetFromMatrix(countData = raw_counts_bck_gs, colData = samplesGsBck2,design = ~originSite_finalSite_experiment) 


# pre-filtering
keep <- rowSums(counts(ddsBck)) >= 10 
ddsBck <- ddsBck[keep,]
keep <- rowSums(counts(ddsBckGs)) >= 10 
ddsBckGs <- ddsBckGs[keep,]
keep <- rowSums(counts(ddsBckGs2)) >= 10 
ddsBckGs2 <- ddsBckGs2[keep,]

# Differential expression analysis
ddsBck<-DESeq(ddsBck)
cbind(resultsNames(ddsBck))
ddsBckGs<-DESeq(ddsBckGs)
cbind(resultsNames(ddsBckGs))
ddsBckGs2<-DESeq(ddsBckGs2)
cbind(resultsNames(ddsBckGs2))


# General 
vsdSept2018Bck = vst(ddsBck,blind=T)
vsdSept2018BckGs = vst(ddsBckGs,blind=T)
vsdSept2018BckGs2 = vst(ddsBckGs2,blind=T)

pcaData = plotPCA(vsdSept2018Bck, intgroup="originSite_finalSite_experiment",returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment)) + 
  geom_point(size = 12) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  geom_point() +
  ggtitle("Home environments - September 2018") +
  theme(text = element_text(size=24),legend.text = element_text(size=20), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 

pcaData = plotPCA(vsdSept2018BckGs, intgroup=c("site","experiment"),returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = site,shape = experiment)) + 
  geom_point(size = 12) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  scale_shape_manual(values = c("circle", "triangle")) +
  geom_point() +
  ggtitle("Home environments & garden same site - September 2018") +
  theme(text = element_text(size=24),legend.text = element_text(size=20), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# Common gm
res_pv_gm_lfc2<-results(ddsBckGs2, contrast=c("originSite_finalSite_experiment","pv_pv_bck","gm_gm_bck"), lfcThreshold = 1,alpha = 0.05)
res_pv_gm_lfc2<-trimFunction(res_pv_gm_lfc2,0.05,1)

res_sp_gm_lfc2<-results(ddsBckGs2, contrast=c("originSite_finalSite_experiment","sp_sp_bck","gm_gm_bck"), lfcThreshold = 1,alpha = 0.05)
res_sp_gm_lfc2<-trimFunction(res_sp_gm_lfc2,0.05,1)

common_gm_lfc2 <- merge(res_pv_gm_lfc2,res_sp_gm_lfc2,by="genes")
common_gm_annot_lfc2 = merge(common_gm_lfc2,A_pfam,by="genes")
common_gm_all_lfc2 <- applyAnnot(common_gm_lfc2,common_gm_annot_lfc2)

newColNames = samplesGsBck2$originSite_finalSite_experiment
vsd_common_gm_sept2018 = heatmapFunction(newColNames,common_gm_lfc2,common_gm_all_lfc2,vsdSept2018BckGs2,"Sept2018 common GM genes")


# Common pv
res_gm_pv_lfc2<-results(ddsBckGs2, contrast=c("originSite_finalSite_experiment","gm_gm_bck","pv_pv_bck"), lfcThreshold = 1,alpha = 0.05)
res_gm_pv_lfc2<-trimFunction(res_gm_pv_lfc2,0.05,1)

res_sp_pv_lfc2<-results(ddsBckGs2, contrast=c("originSite_finalSite_experiment","sp_sp_bck","pv_pv_bck"), lfcThreshold = 1,alpha = 0.05)
res_sp_pv_lfc2<-trimFunction(res_sp_pv_lfc2,0.05,1)

common_pv_lfc2 <- merge(res_gm_pv_lfc2,res_sp_pv_lfc2,by="genes")
common_pv_annot_lfc2 = merge(common_pv_lfc2,A_pfam,by="genes")
common_pv_all_lfc2 <- applyAnnot(common_pv_lfc2,common_pv_annot_lfc2)

newColNames = samplesGsBck2$originSite_finalSite_experiment
vsd_common_pv_sept2018 = heatmapFunction(newColNames,common_pv_lfc2,common_pv_all_lfc2,vsdSept2018BckGs2,"Sept2018 common PV genes")


# Common sp
res_gm_sp_lfc2<-results(ddsBckGs2, contrast=c("originSite_finalSite_experiment","gm_gm_bck","sp_sp_bck"), lfcThreshold = 1,alpha = 0.05)
res_gm_sp_lfc2<-trimFunction(res_gm_sp_lfc2,0.05,1)

res_pv_sp_lfc2<-results(ddsBckGs2, contrast=c("originSite_finalSite_experiment","pv_pv_bck","sp_sp_bck"), lfcThreshold = 1,alpha = 0.05)
res_pv_sp_lfc2<-trimFunction(res_pv_sp_lfc2,0.05,1)

common_sp_lfc2 <- merge(res_gm_sp_lfc2,res_pv_sp_lfc2,by="genes")
common_sp_annot_lfc2 = merge(common_sp_lfc2,A_pfam,by="genes")
common_sp_all_lfc2 <- applyAnnot(common_sp_lfc2,common_sp_annot_lfc2)

newColNames = samplesGsBck2$originSite_finalSite_experiment
vsd_common_sp_sept2018 = heatmapFunction(newColNames,common_sp_lfc2,common_sp_all_lfc2,vsdSept2018BckGs2,"Sept2018 common SP genes")



# Inferences statistics
count_tab_assay <- assay(vsdSept2018Bck)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samplesBck,dist_tab_assay ~ originSite_finalSite_experiment, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesBck$originSite_finalSite_experiment)
pair.mod
mod <- betadisper(dist_tab_assay,samplesBck$originSite_finalSite_experiment)
mod.HSD <- TukeyHSD(mod)
mod.HSD


count_tab_assay <- assay(vsdSept2018BckGs2)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samplesGsBck2,dist_tab_assay ~ originSite_finalSite_experiment, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesGsBck2$originSite_finalSite_experiment)
pair.mod
mod <- betadisper(dist_tab_assay,samplesGsBck2$originSite_finalSite_experiment)
mod.HSD <- TukeyHSD(mod)
mod.HSD


count_tab_assay <- assay(vsdSept2018BckGs)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samplesGsBck,dist_tab_assay ~ site + experiment, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesGsBck$site)
pair.mod
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesGsBck$experiment)
pair.mod
mod <- betadisper(dist_tab_assay,samplesGsBck$site)
mod.HSD <- TukeyHSD(mod)
mod.HSD


# May2018 & Sept2018 - Background
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
A_pfam<-read.csv('Astroides_pfam.csv',header=T,sep=',')
samplesBck<-read.table('tximport_design_background2018.txt',header=T)
dataPathBck<-'/Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/bck'
setwd(dataPathBck)
data<-list.files(pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts.files <- lapply(data, read.table, skip = 4)
raw_counts <- as.data.frame(sapply(counts.files, function(x) x[ , 2]))
data <- gsub( "Users/mmeynadier/Documents/PhD/species/Astroides/analysis/STARmapping/teixido/adult/bck", "", data )
data <- gsub( "_ReadsPerGene.out.tab", "", data )
data <- gsub( "./", "", data )
colnames(raw_counts) <- data
row.names(raw_counts) <- counts.files[[1]]$V1
samplesOrder = c(4,5,6,13,14,15,1,2,3,10,11,12,7,8,9,16,17,18)
raw_counts <- raw_counts[,samplesOrder]

# DDS object
dds<-DESeqDataSetFromMatrix(countData = raw_counts, colData = samplesBck,design = ~originSite_finalSite_experiment + dataset)

# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]

# Differential expression analysis
dds<-DESeq(dds)
cbind(resultsNames(dds))

# General 
vsd = vst(dds,blind=T)

pcaData = plotPCA(vsd, intgroup=c("originSite_finalSite_experiment","dataset"),returnData=TRUE,ntop=200)
percentVar = round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, colour = originSite_finalSite_experiment,shape=dataset)) + 
  geom_point(size = 12) + theme_bw() + 
  scale_color_manual(values = c("#ff4040", "#00008B","#6495ED")) +
  scale_shape_manual(values = c("circle", "triangle")) +
  geom_point() +
  ggtitle("Home environments - May 2018 & September 2018") +
  theme(text = element_text(size=24),legend.text = element_text(size=20), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 


# Common gm
res_pv_gm_lfc2<-results(dds, contrast=c("originSite_finalSite_experiment","pv_pv_bck","gm_gm_bck"), lfcThreshold = 1,alpha = 0.05)
res_pv_gm_lfc2<-trimFunction(res_pv_gm_lfc2,0.05,1)

res_sp_gm_lfc2<-results(dds, contrast=c("originSite_finalSite_experiment","sp_sp_bck","gm_gm_bck"), lfcThreshold = 1,alpha = 0.05)
res_sp_gm_lfc2<-trimFunction(res_sp_gm_lfc2,0.05,1)

common_gm_lfc2 <- merge(res_pv_gm_lfc2,res_sp_gm_lfc2,by="genes")
common_gm_annot_lfc2 = merge(common_gm_lfc2,A_pfam,by="genes")
common_gm_all_lfc2 <- applyAnnot(common_gm_lfc2,common_gm_annot_lfc2)

newColNames = samplesBck$originSite_finalSite_experiment
vsd_common_gm = heatmapFunction(newColNames,common_gm_lfc2,common_gm_all_lfc2,vsd,"May2018 & Sept2018 common GM genes")


# Common pv
res_gm_pv_lfc2<-results(dds, contrast=c("originSite_finalSite_experiment","gm_gm_bck","pv_pv_bck"), lfcThreshold = 1,alpha = 0.05)
res_gm_pv_lfc2<-trimFunction(res_gm_pv_lfc2,0.05,1)

res_sp_pv_lfc2<-results(dds, contrast=c("originSite_finalSite_experiment","sp_sp_bck","pv_pv_bck"), lfcThreshold = 1,alpha = 0.05)
res_sp_pv_lfc2<-trimFunction(res_sp_pv_lfc2,0.05,1)

common_pv_lfc2 <- merge(res_gm_pv_lfc2,res_sp_pv_lfc2,by="genes")
common_pv_annot_lfc2 = merge(common_pv_lfc2,A_pfam,by="genes")
common_pv_all_lfc2 <- applyAnnot(common_pv_lfc2,common_pv_annot_lfc2)

newColNames = samplesBck$originSite_finalSite_experiment
vsd_common_pv_ = heatmapFunction(newColNames,common_pv_lfc2,common_pv_all_lfc2,vsd,"May2018 & Sept2018 common PV genes")


# Common sp
res_gm_sp_lfc2<-results(dds, contrast=c("originSite_finalSite_experiment","gm_gm_bck","sp_sp_bck"), lfcThreshold = 1,alpha = 0.05)
res_gm_sp_lfc2<-trimFunction(res_gm_sp_lfc2,0.05,1)

res_pv_sp_lfc2<-results(dds, contrast=c("originSite_finalSite_experiment","pv_pv_bck","sp_sp_bck"), lfcThreshold = 1,alpha = 0.05)
res_pv_sp_lfc2<-trimFunction(res_pv_sp_lfc2,0.05,1)

common_sp_lfc2 <- merge(res_gm_sp_lfc2,res_pv_sp_lfc2,by="genes")
common_sp_annot_lfc2 = merge(common_sp_lfc2,A_pfam,by="genes")
common_sp_all_lfc2 <- applyAnnot(common_sp_lfc2,common_sp_annot_lfc2)

newColNames = samplesBck$originSite_finalSite_experiment
vsd_common_sp = heatmapFunction(newColNames,common_sp_lfc2,common_sp_all_lfc2,vsd,"May2018 & Sept2018 common SP genes")


# Inferences statistics
count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
ad<-adonis2(data=samplesBck,dist_tab_assay ~ originSite_finalSite_experiment + dataset, method="euclidian")
ad
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesBck$originSite_finalSite_experiment)
pair.mod
pair.mod<-pairwise.adonis(dist_tab_assay,factors=samplesBck$dataset)
pair.mod
mod <- betadisper(dist_tab_assay,samplesBck$originSite_finalSite_experiment)
mod.HSD <- TukeyHSD(mod)
mod.HSD
mod <- betadisper(dist_tab_assay,samplesBck$dataset)
mod.HSD <- TukeyHSD(mod)
mod.HSD
