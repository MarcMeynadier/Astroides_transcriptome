DE_Astroides_juvenile
================
Marc Meynadier
4/11/2022

``` r
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
```

    ## Le chargement a nécessité le package : DESeq2

    ## Le chargement a nécessité le package : S4Vectors

    ## Le chargement a nécessité le package : stats4

    ## Le chargement a nécessité le package : BiocGenerics

    ## 
    ## Attachement du package : 'BiocGenerics'

    ## Les objets suivants sont masqués depuis 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## Les objets suivants sont masqués depuis 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## 
    ## Attachement du package : 'S4Vectors'

    ## Les objets suivants sont masqués depuis 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Le chargement a nécessité le package : IRanges

    ## Le chargement a nécessité le package : GenomicRanges

    ## Le chargement a nécessité le package : GenomeInfoDb

    ## Le chargement a nécessité le package : SummarizedExperiment

    ## Le chargement a nécessité le package : MatrixGenerics

    ## Le chargement a nécessité le package : matrixStats

    ## 
    ## Attachement du package : 'MatrixGenerics'

    ## Les objets suivants sont masqués depuis 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Le chargement a nécessité le package : Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attachement du package : 'Biobase'

    ## L'objet suivant est masqué depuis 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## Les objets suivants sont masqués depuis 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## Le chargement a nécessité le package : devtools

    ## Le chargement a nécessité le package : usethis

    ## Le chargement a nécessité le package : BiocManager

    ## Bioconductor version '3.14' is out-of-date; the current release version '3.15'
    ##   is available with R version '4.2'; see https://bioconductor.org/install

    ## 
    ## Attachement du package : 'BiocManager'

    ## L'objet suivant est masqué depuis 'package:devtools':
    ## 
    ##     install

    ## Le chargement a nécessité le package : ggplot2

    ## Le chargement a nécessité le package : ggrepel

    ## Le chargement a nécessité le package : markdown

    ## Le chargement a nécessité le package : RColorBrewer

    ## Le chargement a nécessité le package : genefilter

    ## 
    ## Attachement du package : 'genefilter'

    ## Les objets suivants sont masqués depuis 'package:MatrixGenerics':
    ## 
    ##     rowSds, rowVars

    ## Les objets suivants sont masqués depuis 'package:matrixStats':
    ## 
    ##     rowSds, rowVars

    ## Le chargement a nécessité le package : gplots

    ## 
    ## Attachement du package : 'gplots'

    ## L'objet suivant est masqué depuis 'package:IRanges':
    ## 
    ##     space

    ## L'objet suivant est masqué depuis 'package:S4Vectors':
    ## 
    ##     space

    ## L'objet suivant est masqué depuis 'package:stats':
    ## 
    ##     lowess

    ## Le chargement a nécessité le package : vegan

    ## Le chargement a nécessité le package : permute

    ## 
    ## Attachement du package : 'permute'

    ## L'objet suivant est masqué depuis 'package:devtools':
    ## 
    ##     check

    ## Le chargement a nécessité le package : lattice

    ## This is vegan 2.5-7

    ## Le chargement a nécessité le package : dplyr

    ## 
    ## Attachement du package : 'dplyr'

    ## L'objet suivant est masqué depuis 'package:Biobase':
    ## 
    ##     combine

    ## L'objet suivant est masqué depuis 'package:matrixStats':
    ## 
    ##     count

    ## Les objets suivants sont masqués depuis 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## L'objet suivant est masqué depuis 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## Les objets suivants sont masqués depuis 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## Les objets suivants sont masqués depuis 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## Les objets suivants sont masqués depuis 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## Les objets suivants sont masqués depuis 'package:stats':
    ## 
    ##     filter, lag

    ## Les objets suivants sont masqués depuis 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
#BiocManager::install('tximport', force = TRUE)
#BiocManager::install('apeglm')
#BiocManager::install('ashr')
#BiocManager::install("EnhancedVolcano")
#BiocManager::install('limma')
#devtools::install_github('cran/GMD')
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
```

    ## Skipping install of 'ggvenn' from a github remote, the SHA1 (b7ff54ba) has not changed since last install.
    ##   Use `force = TRUE` to force installation

``` r
library('ggvenn')
```

    ## Le chargement a nécessité le package : grid

``` r
library('tximport')
library('apeglm')
library('ashr')
library('EnhancedVolcano')
```

    ## Registered S3 methods overwritten by 'ggalt':
    ##   method                  from   
    ##   grid.draw.absoluteGrob  ggplot2
    ##   grobHeight.absoluteGrob ggplot2
    ##   grobWidth.absoluteGrob  ggplot2
    ##   grobX.absoluteGrob      ggplot2
    ##   grobY.absoluteGrob      ggplot2

``` r
library('limma')
```

    ## 
    ## Attachement du package : 'limma'

    ## L'objet suivant est masqué depuis 'package:DESeq2':
    ## 
    ##     plotMA

    ## L'objet suivant est masqué depuis 'package:BiocGenerics':
    ## 
    ##     plotMA

``` r
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
```

    ## ℹ SHA-1 hash of file is 015fc0457e61e3e93a903e69a24d96d2dac7b9fb

``` r
# Working environment 
scriptPath<-dirname(rstudioapi::getSourceEditorContext()$path)
setwd(scriptPath)
samples<-read.table('tximport_design_juvenile.txt',header=T)
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
filesNatSim<-paste0(samplesNatSim$samples,'.tsv')
names(files)<-samples$samples
names(filesNatSim)<-samplesNatSim$samples
txi<-tximport(files = files,type='kallisto',tx2gene = tx2gene)
```

    ## Note: importing `abundance.h5` is typically faster than `abundance.tsv`

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 
    ## transcripts missing from tx2gene: 1
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
txiNatSim<-tximport(files = filesNatSim,type='kallisto',tx2gene = tx2gene)
```

    ## Note: importing `abundance.h5` is typically faster than `abundance.tsv`
    ## reading in files with read_tsv
    ## 1 2 3 4 5 6 7 8 
    ## transcripts missing from tx2gene: 1
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
names(txi)
```

    ## [1] "abundance"           "counts"              "length"             
    ## [4] "countsFromAbundance"

``` r
names(txiNatSim)
```

    ## [1] "abundance"           "counts"              "length"             
    ## [4] "countsFromAbundance"

``` r
head(txi$counts)
```

    ##                    abundance_juvenile_rep1_gm_amb_R1CA1_
    ## TRINITY_DN0_c0_g1                                725.000
    ## TRINITY_DN0_c0_g2                                252.783
    ## TRINITY_DN0_c1_g1                                148.439
    ## TRINITY_DN1_c0_g1                               5981.350
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1721.830
    ##                    abundance_juvenile_rep2_gm_amb_R2DA8_
    ## TRINITY_DN0_c0_g1                               1403.000
    ## TRINITY_DN0_c0_g2                                248.403
    ## TRINITY_DN0_c1_g1                                104.596
    ## TRINITY_DN1_c0_g1                               6103.210
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1704.760
    ##                    abundance_juvenile_rep3_gm_amb_R3AA9_
    ## TRINITY_DN0_c0_g1                              1247.0000
    ## TRINITY_DN0_c0_g2                               140.4600
    ## TRINITY_DN0_c1_g1                                64.3743
    ## TRINITY_DN1_c0_g1                              7743.7300
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1980.4300
    ##                    abundance_juvenile_rep4_gm_amb_R4BA4_
    ## TRINITY_DN0_c0_g1                               1390.000
    ## TRINITY_DN0_c0_g2                                267.591
    ## TRINITY_DN0_c1_g1                                143.146
    ## TRINITY_DN1_c0_g1                               7346.000
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2026.520
    ##                    abundance_juvenile_rep1_gm_low_R1BA8_
    ## TRINITY_DN0_c0_g1                               637.0000
    ## TRINITY_DN0_c0_g2                               137.4180
    ## TRINITY_DN0_c1_g1                                75.0977
    ## TRINITY_DN1_c0_g1                              8165.0000
    ## TRINITY_DN1_c1_g1                                 0.0000
    ## TRINITY_DN10_c0_g1                             1309.2800
    ##                    abundance_juvenile_rep2_gm_low_R2EA3_
    ## TRINITY_DN0_c0_g1                                799.000
    ## TRINITY_DN0_c0_g2                                338.944
    ## TRINITY_DN0_c1_g1                                108.556
    ## TRINITY_DN1_c0_g1                               9692.170
    ## TRINITY_DN1_c1_g1                                  4.000
    ## TRINITY_DN10_c0_g1                              1149.120
    ##                    abundance_juvenile_rep3_gm_low_R3DA1_
    ## TRINITY_DN0_c0_g1                               1099.000
    ## TRINITY_DN0_c0_g2                                332.781
    ## TRINITY_DN0_c1_g1                                145.105
    ## TRINITY_DN1_c0_g1                               6556.000
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1808.160
    ##                    abundance_juvenile_rep4_gm_low_R4EA1_
    ## TRINITY_DN0_c0_g1                               1205.000
    ## TRINITY_DN0_c0_g2                                210.875
    ## TRINITY_DN0_c1_g1                                131.151
    ## TRINITY_DN1_c0_g1                               7167.150
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1575.450
    ##                    abundance_juvenile_rep1_gm_ext_R1EA7_
    ## TRINITY_DN0_c0_g1                                862.000
    ## TRINITY_DN0_c0_g2                                249.846
    ## TRINITY_DN0_c1_g1                                180.206
    ## TRINITY_DN1_c0_g1                               8382.000
    ## TRINITY_DN1_c1_g1                                  2.000
    ## TRINITY_DN10_c0_g1                              1785.440
    ##                    abundance_juvenile_rep2_gm_ext_R2AA2_
    ## TRINITY_DN0_c0_g1                               792.0000
    ## TRINITY_DN0_c0_g2                                98.9848
    ## TRINITY_DN0_c1_g1                               133.4120
    ## TRINITY_DN1_c0_g1                              4871.0000
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1624.2100
    ##                    abundance_juvenile_rep3_gm_ext_R3CA1_
    ## TRINITY_DN0_c0_g1                               1040.000
    ## TRINITY_DN0_c0_g2                                160.000
    ## TRINITY_DN0_c1_g1                                139.177
    ## TRINITY_DN1_c0_g1                               7980.000
    ## TRINITY_DN1_c1_g1                                  2.000
    ## TRINITY_DN10_c0_g1                              1813.350
    ##                    abundance_juvenile_rep4_gm_ext_R4BA15_
    ## TRINITY_DN0_c0_g1                                 707.000
    ## TRINITY_DN0_c0_g2                                 319.870
    ## TRINITY_DN0_c1_g1                                 123.098
    ## TRINITY_DN1_c0_g1                                8087.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               1589.580
    ##                    abundance_juvenile_rep1_sp_amb_R1AA8_
    ## TRINITY_DN0_c0_g1                               1192.000
    ## TRINITY_DN0_c0_g2                                117.032
    ## TRINITY_DN0_c1_g1                                102.447
    ## TRINITY_DN1_c0_g1                               7894.600
    ## TRINITY_DN1_c1_g1                                  6.000
    ## TRINITY_DN10_c0_g1                              1612.620
    ##                    abundance_juvenile_rep2_sp_amb_R2BA7_
    ## TRINITY_DN0_c0_g1                               947.0000
    ## TRINITY_DN0_c0_g2                               224.4870
    ## TRINITY_DN0_c1_g1                                97.6384
    ## TRINITY_DN1_c0_g1                              6051.0000
    ## TRINITY_DN1_c1_g1                                 2.0000
    ## TRINITY_DN10_c0_g1                             1869.9000
    ##                    abundance_juvenile_rep3_sp_amb_R3EA10_
    ## TRINITY_DN0_c0_g1                                1278.000
    ## TRINITY_DN0_c0_g2                                 224.125
    ## TRINITY_DN0_c1_g1                                 120.269
    ## TRINITY_DN1_c0_g1                                6596.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               2252.890
    ##                    abundance_juvenile_rep4_sp_amb_R4AA12_
    ## TRINITY_DN0_c0_g1                               1054.0000
    ## TRINITY_DN0_c0_g2                                168.0360
    ## TRINITY_DN0_c1_g1                                 97.0586
    ## TRINITY_DN1_c0_g1                               7395.0000
    ## TRINITY_DN1_c1_g1                                  0.0000
    ## TRINITY_DN10_c0_g1                              2042.4100
    ##                    abundance_juvenile_rep1_sp_low_R1DA1_
    ## TRINITY_DN0_c0_g1                               923.0000
    ## TRINITY_DN0_c0_g2                               222.0880
    ## TRINITY_DN0_c1_g1                                90.0538
    ## TRINITY_DN1_c0_g1                              6924.0000
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1762.3000
    ##                    abundance_juvenile_rep2_sp_low_R2CA9_
    ## TRINITY_DN0_c0_g1                                924.000
    ## TRINITY_DN0_c0_g2                                208.798
    ## TRINITY_DN0_c1_g1                                132.564
    ## TRINITY_DN1_c0_g1                               6535.180
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1409.980
    ##                    abundance_juvenile_rep3_sp_low_R3BA6_
    ## TRINITY_DN0_c0_g1                               1397.000
    ## TRINITY_DN0_c0_g2                                230.733
    ## TRINITY_DN0_c1_g1                                124.001
    ## TRINITY_DN1_c0_g1                               6544.980
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2067.810
    ##                    abundance_juvenile_rep4_sp_low_R4CA7_
    ## TRINITY_DN0_c0_g1                               1041.000
    ## TRINITY_DN0_c0_g2                                 93.000
    ## TRINITY_DN0_c1_g1                                105.654
    ## TRINITY_DN1_c0_g1                               9480.650
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1832.770

``` r
head(txiNatSim$counts)
```

    ##                    abundance_juvenile_rep1_gm_low_R1BA8_
    ## TRINITY_DN0_c0_g1                               637.0000
    ## TRINITY_DN0_c0_g2                               137.4180
    ## TRINITY_DN0_c1_g1                                75.0977
    ## TRINITY_DN1_c0_g1                              8165.0000
    ## TRINITY_DN1_c1_g1                                 0.0000
    ## TRINITY_DN10_c0_g1                             1309.2800
    ##                    abundance_juvenile_rep2_gm_low_R2EA3_
    ## TRINITY_DN0_c0_g1                                799.000
    ## TRINITY_DN0_c0_g2                                338.944
    ## TRINITY_DN0_c1_g1                                108.556
    ## TRINITY_DN1_c0_g1                               9692.170
    ## TRINITY_DN1_c1_g1                                  4.000
    ## TRINITY_DN10_c0_g1                              1149.120
    ##                    abundance_juvenile_rep3_gm_low_R3DA1_
    ## TRINITY_DN0_c0_g1                               1099.000
    ## TRINITY_DN0_c0_g2                                332.781
    ## TRINITY_DN0_c1_g1                                145.105
    ## TRINITY_DN1_c0_g1                               6556.000
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1808.160
    ##                    abundance_juvenile_rep4_gm_low_R4EA1_
    ## TRINITY_DN0_c0_g1                               1205.000
    ## TRINITY_DN0_c0_g2                                210.875
    ## TRINITY_DN0_c1_g1                                131.151
    ## TRINITY_DN1_c0_g1                               7167.150
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1575.450
    ##                    abundance_juvenile_rep1_sp_amb_R1AA8_
    ## TRINITY_DN0_c0_g1                               1192.000
    ## TRINITY_DN0_c0_g2                                117.032
    ## TRINITY_DN0_c1_g1                                102.447
    ## TRINITY_DN1_c0_g1                               7894.600
    ## TRINITY_DN1_c1_g1                                  6.000
    ## TRINITY_DN10_c0_g1                              1612.620
    ##                    abundance_juvenile_rep2_sp_amb_R2BA7_
    ## TRINITY_DN0_c0_g1                               947.0000
    ## TRINITY_DN0_c0_g2                               224.4870
    ## TRINITY_DN0_c1_g1                                97.6384
    ## TRINITY_DN1_c0_g1                              6051.0000
    ## TRINITY_DN1_c1_g1                                 2.0000
    ## TRINITY_DN10_c0_g1                             1869.9000
    ##                    abundance_juvenile_rep3_sp_amb_R3EA10_
    ## TRINITY_DN0_c0_g1                                1278.000
    ## TRINITY_DN0_c0_g2                                 224.125
    ## TRINITY_DN0_c1_g1                                 120.269
    ## TRINITY_DN1_c0_g1                                6596.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               2252.890
    ##                    abundance_juvenile_rep4_sp_amb_R4AA12_
    ## TRINITY_DN0_c0_g1                               1054.0000
    ## TRINITY_DN0_c0_g2                                168.0360
    ## TRINITY_DN0_c1_g1                                 97.0586
    ## TRINITY_DN1_c0_g1                               7395.0000
    ## TRINITY_DN1_c1_g1                                  0.0000
    ## TRINITY_DN10_c0_g1                              2042.4100

``` r
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site+pH)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

``` r
ddsNatSim<-DESeqDataSetFromTximport(txiNatSim,colData=samplesNatSim,design= ~site_pH)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

``` r
# pre-filtering
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,]
keep <- rowSums(counts(ddsNatSim)) >= 10 
ddsNatSim <- ddsNatSim[keep,]

# Differential expression analysis
dds<-DESeq(dds)
```

    ## estimating size factors
    ## using 'avgTxLength' from assays(dds), correcting for library size
    ## estimating dispersions
    ## gene-wise dispersion estimates
    ## mean-dispersion relationship
    ## final dispersion estimates
    ## fitting model and testing

``` r
ddsNatSim<-DESeq(ddsNatSim)
```

    ## estimating size factors
    ## using 'avgTxLength' from assays(dds), correcting for library size
    ## estimating dispersions
    ## gene-wise dispersion estimates
    ## mean-dispersion relationship
    ## final dispersion estimates
    ## fitting model and testing

``` r
cbind(resultsNames(dds))
```

    ##      [,1]           
    ## [1,] "Intercept"    
    ## [2,] "site_sp_vs_gm"
    ## [3,] "pH_ext_vs_amb"
    ## [4,] "pH_low_vs_amb"

``` r
cbind(resultsNames(ddsNatSim))
```

    ##      [,1]                      
    ## [1,] "Intercept"               
    ## [2,] "site_pH_sp_amb_vs_gm_low"

``` r
sp_VS_gm<-results(dds, contrast=c("site","sp","gm"), alpha = 0.05)
amb_VS_ext<-results(dds, contrast=c("pH","amb","ext"), alpha = 0.05)
amb_VS_low<-results(dds, contrast=c("pH","amb","low"), alpha = 0.05)
low_VS_ext<-results(dds, contrast=c("pH","low","ext"), alpha = 0.05)
sp_amb_VS_gm_low_natSim<-results(ddsNatSim,contrast=c("site_pH","sp_amb","gm_low"),alpha = 0.05)
summary(sp_VS_gm)
```

    ## 
    ## out of 93585 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 378, 0.4%
    ## LFC < 0 (down)     : 509, 0.54%
    ## outliers [1]       : 583, 0.62%
    ## low counts [2]     : 7227, 7.7%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(amb_VS_ext)
```

    ## 
    ## out of 93585 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 31, 0.033%
    ## LFC < 0 (down)     : 12, 0.013%
    ## outliers [1]       : 583, 0.62%
    ## low counts [2]     : 1813, 1.9%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(amb_VS_low)
```

    ## 
    ## out of 93585 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 6, 0.0064%
    ## LFC < 0 (down)     : 5, 0.0053%
    ## outliers [1]       : 583, 0.62%
    ## low counts [2]     : 21617, 23%
    ## (mean count < 5)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(low_VS_ext)
```

    ## 
    ## out of 93585 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 38, 0.041%
    ## LFC < 0 (down)     : 11, 0.012%
    ## outliers [1]       : 583, 0.62%
    ## low counts [2]     : 18018, 19%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(sp_amb_VS_gm_low_natSim)
```

    ## 
    ## out of 87637 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 355, 0.41%
    ## LFC < 0 (down)     : 371, 0.42%
    ## outliers [1]       : 680, 0.78%
    ## low counts [2]     : 11893, 14%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
# Exploring the results

# Results sp VS gm

#MA-plot
DESeq2::plotMA(sp_VS_gm,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nsp_VS_gm")
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# Volcano plot
pCutoff = 0.05
FCcutoff = 1.0
EnhancedVolcano(data.frame(sp_VS_gm), lab = rownames(data.frame(sp_VS_gm)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp and gm",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_VS_gm), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
# Results ext VS amb

#MA-plot
DESeq2::plotMA(amb_VS_ext,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\namb_VS_ext")
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
# Volcano plot
EnhancedVolcano(data.frame(amb_VS_ext), lab = rownames(data.frame(amb_VS_ext)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between ext and amb",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(amb_VS_ext), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->

``` r
# Results low VS amb

#MA-plot
DESeq2::plotMA(amb_VS_low,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\namb_VS_low")
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->

``` r
# Volcano plot
EnhancedVolcano(data.frame(amb_VS_low), lab = rownames(data.frame(amb_VS_low)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between low and amb",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(amb_VS_low), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->

``` r
# Results low VS ext

#MA-plot
DESeq2::plotMA(low_VS_ext,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nlow_VS_ext")
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-7.png)<!-- -->

``` r
# Volcano plot
EnhancedVolcano(data.frame(low_VS_ext), lab = rownames(data.frame(low_VS_ext)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between low and ext",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(low_VS_ext), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-8.png)<!-- -->

``` r
# Results natural simulation

#MA-plot
DESeq2::plotMA(sp_amb_VS_gm_low_natSim,ylim=c(-50,50),main="MA-plot for the shrunken log2 fold changes\nsp_amb_VS_gm_low")
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-9.png)<!-- -->

``` r
# Volcano plot
EnhancedVolcano(data.frame(sp_amb_VS_gm_low_natSim), lab = rownames(data.frame(sp_amb_VS_gm_low_natSim)), x = 'log2FoldChange', y = 'padj',
                xlab = bquote(~Log[2]~ 'fold change'), ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = pCutoff, FCcutoff = FCcutoff, pointSize = 1.0, labSize = 2.0,
                title = "Volcano plot", subtitle = "Contrast between sp_amb and gm_low",
                caption = paste0('log2 FC cutoff: ', FCcutoff, '; p-value cutoff: ', pCutoff, '\nTotal = ', nrow(sp_amb_VS_gm_low_natSim), ' variables'),
                legendLabels=c('NS','Log2 FC','Adjusted p-value', 'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom', legendLabSize = 14, legendIconSize = 5.0)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-10.png)<!-- -->

``` r
# Principal Component Analysis
vsd = vst(dds,blind=T)

pcaData = plotPCA(vsd, intgroup=c("site","pH"), 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, colour = site, shape = pH)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040","#6495ED")) +
  scale_shape_manual(values = c("triangle","circle","square")) +
  geom_point() +
  ggtitle("Principal Component Analysis of juvenile corals", subtitle = "Juvenile dataset") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-11.png)<!-- -->

``` r
vsdNatSim = vst(ddsNatSim,blind=T)

pcaData = plotPCA(vsdNatSim, intgroup="site_pH", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, colour = site_pH)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#ff4040","#6495ED")) +
  geom_point() +
  ggtitle("Principal Component Analysis of juvenile corals", subtitle = "Juvenile dataset - Natural conditions simulation") +
  theme(text = element_text(size=14),legend.text = element_text(size=12), legend.position = 'bottom') +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) 
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-12.png)<!-- -->

``` r
# Venn diagramm 
resOrdered_sp_VS_gm <- sp_VS_gm[order(sp_VS_gm$padj),]
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

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#009E73"),
  stroke_size = 0.5, set_name_size = 4
)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-13.png)<!-- -->

``` r
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
heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",key.title = "none",
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7,density.info="none",
          xlab="sampling sites",ylab="proteins associated to genes",Colv=NA,margins = c(5, 7))
```

    ## Warning in heatmap.2(assay(vsdCandidate)[topVarGenesVsd, ], trace = "none", :
    ## Discrepancy: Colv is FALSE, while dendrogram is `both'. Omitting column
    ## dendogram.

``` r
main='Differential expression of 50 most expressed candidates genes\n\nJuveniles'
title(main, cex.main = 0.7)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-14.png)<!-- -->

``` r
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
heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",key.title = "none",
          col=colorRampPalette(rev(brewer.pal(11,"PuOr")))(255), cexRow=0.6, cexCol=0.7,density.info="none",
          xlab="sampling sites",ylab="proteins associated to genes",Colv=NA,margins = c(5, 7))
```

    ## Warning in heatmap.2(assay(vsdCandidate)[topVarGenesVsd, ], trace = "none", :
    ## Discrepancy: Colv is FALSE, while dendrogram is `both'. Omitting column
    ## dendogram.

``` r
main='Differential expression of 50 most expressed candidates genes\n\nJuveniles - Natural conditions simulation'
title(main, cex.main = 0.7)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-15.png)<!-- -->

``` r
# Inferences statistics

count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site + pH, method="euclidian")
```

    ## 
    ## Call:
    ## adonis(formula = dist_tab_assay ~ site + pH, data = samples,      method = "euclidian") 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## site       1     36449   36449  2.0478 0.10260  0.001 ***
    ## pH         2     34016   17008  0.9556 0.09576  0.676    
    ## Residuals 16    284777   17799         0.80164           
    ## Total     19    355242                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(betadisper(dist_tab_assay,samples$site))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq F value Pr(>F)
    ## Groups     1  343.43  343.43    2.98 0.1014
    ## Residuals 18 2074.40  115.24

``` r
anova(betadisper(dist_tab_assay,samples$pH))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## Groups     2 1088.5  544.23  6.2804 0.009073 **
    ## Residuals 17 1473.2   86.66                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
count_tab_assay <- assay(vsdNatSim)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samplesNatSim,dist_tab_assay ~ site_pH, method="euclidian")
```

    ## 
    ## Call:
    ## adonis(formula = dist_tab_assay ~ site_pH, data = samplesNatSim,      method = "euclidian") 
    ## 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
    ## site_pH    1     26258   26258  1.5539 0.20571  0.042 *
    ## Residuals  6    101388   16898         0.79429         
    ## Total      7    127646                 1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(betadisper(dist_tab_assay,samplesNatSim$site_pH))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## Groups     1 288.87 288.868  3.3437 0.1172
    ## Residuals  6 518.35  86.392

``` r
# Exporting results
resOrdered_sp_amb_VS_gm_low_natSim <- sp_amb_VS_gm_low_natSim[order(sp_amb_VS_gm_low_natSim$padj),]
resOrderedDF_sp_amb_VS_gm_low_natSim <- as.data.frame(resOrdered_sp_amb_VS_gm_low_natSim)

write.csv(resOrderedDF_sp_VS_gm, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_sp_VS_gm.csv',sep=''))
write.csv(resOrderedDF_amb_VS_ext, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_amb_VS_ext.csv',sep=''))
write.csv(resOrderedDF_amb_VS_low, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_amb_VS_low.csv',sep=''))
write.csv(resOrderedDF_low_VS_ext, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_low_VS_ext.csv',sep=''))
write.csv(resOrderedDF_sp_amb_VS_gm_low_natSim, file = paste(scriptPath,'/data/net/7_deseq2/adultTranscriptome/juvenile/DESeq2_results_juvenile_sp_amb_VS_gm_low_natSim.csv',sep=''))

sessionInfo()
```

    ## R version 4.1.2 (2021-11-01)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] fr_FR.UTF-8/fr_FR.UTF-8/fr_FR.UTF-8/C/fr_FR.UTF-8/fr_FR.UTF-8
    ## 
    ## attached base packages:
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] limma_3.50.1                EnhancedVolcano_1.12.0     
    ##  [3] ashr_2.2-54                 apeglm_1.16.0              
    ##  [5] tximport_1.22.0             ggvenn_0.1.9               
    ##  [7] dplyr_1.0.8                 vegan_2.5-7                
    ##  [9] lattice_0.20-45             permute_0.9-7              
    ## [11] gplots_3.1.1                genefilter_1.76.0          
    ## [13] RColorBrewer_1.1-3          markdown_1.1               
    ## [15] ggrepel_0.9.1               ggplot2_3.3.5              
    ## [17] BiocManager_1.30.16         devtools_2.4.3             
    ## [19] usethis_2.1.5               DESeq2_1.34.0              
    ## [21] SummarizedExperiment_1.24.0 Biobase_2.54.0             
    ## [23] MatrixGenerics_1.6.0        matrixStats_0.61.0         
    ## [25] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
    ## [27] IRanges_2.28.0              S4Vectors_0.32.3           
    ## [29] BiocGenerics_0.40.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] plyr_1.8.6             splines_4.1.2          BiocParallel_1.28.3   
    ##   [4] digest_0.6.29          invgamma_1.1           htmltools_0.5.2       
    ##   [7] SQUAREM_2021.1         fansi_1.0.3            magrittr_2.0.3        
    ##  [10] memoise_2.0.1          cluster_2.1.2          tzdb_0.2.0            
    ##  [13] remotes_2.4.2          Biostrings_2.62.0      readr_2.1.2           
    ##  [16] annotate_1.72.0        extrafont_0.17         vroom_1.5.7           
    ##  [19] extrafontdb_1.0        bdsmatrix_1.3-4        prettyunits_1.1.1     
    ##  [22] colorspace_2.0-3       blob_1.2.2             xfun_0.30             
    ##  [25] callr_3.7.0            crayon_1.5.1           RCurl_1.98-1.6        
    ##  [28] survival_3.3-1         glue_1.6.2             gtable_0.3.0          
    ##  [31] zlibbioc_1.40.0        XVector_0.34.0         DelayedArray_0.20.0   
    ##  [34] proj4_1.0-11           pkgbuild_1.3.1         Rttf2pt1_1.3.10       
    ##  [37] maps_3.4.0             scales_1.2.0           mvtnorm_1.1-3         
    ##  [40] DBI_1.1.2              Rcpp_1.0.8.2           xtable_1.8-4          
    ##  [43] emdbook_1.3.12         bit_4.0.4              truncnorm_1.0-8       
    ##  [46] httr_1.4.2             ellipsis_0.3.2         farver_2.1.0          
    ##  [49] pkgconfig_2.0.3        XML_3.99-0.9           locfit_1.5-9.5        
    ##  [52] utf8_1.2.2             labeling_0.4.2         tidyselect_1.1.2      
    ##  [55] rlang_1.0.2            AnnotationDbi_1.56.2   munsell_0.5.0         
    ##  [58] tools_4.1.2            cachem_1.0.6           cli_3.2.0             
    ##  [61] generics_0.1.2         RSQLite_2.2.10         evaluate_0.15         
    ##  [64] stringr_1.4.0          fastmap_1.1.0          yaml_2.3.5            
    ##  [67] processx_3.5.2         knitr_1.37             bit64_4.0.5           
    ##  [70] fs_1.5.2               caTools_1.18.2         purrr_0.3.4           
    ##  [73] KEGGREST_1.34.0        nlme_3.1-155           ash_1.0-15            
    ##  [76] ggrastr_1.0.1          brio_1.1.3             compiler_4.1.2        
    ##  [79] rstudioapi_0.13        beeswarm_0.4.0         curl_4.3.2            
    ##  [82] png_0.1-7              testthat_3.1.2         tibble_3.1.6          
    ##  [85] geneplotter_1.72.0     stringi_1.7.6          ps_1.6.0              
    ##  [88] desc_1.4.1             ggalt_0.4.0            Matrix_1.4-0          
    ##  [91] vctrs_0.4.1            pillar_1.7.0           lifecycle_1.0.1       
    ##  [94] bitops_1.0-7           irlba_2.3.5            R6_2.5.1              
    ##  [97] KernSmooth_2.23-20     vipor_0.4.5            sessioninfo_1.2.2     
    ## [100] MASS_7.3-55            gtools_3.9.2           assertthat_0.2.1      
    ## [103] pkgload_1.2.4          rprojroot_2.0.2        withr_2.5.0           
    ## [106] GenomeInfoDbData_1.2.7 mgcv_1.8-39            parallel_4.1.2        
    ## [109] hms_1.1.1              coda_0.19-4            rmarkdown_2.13        
    ## [112] mixsqp_0.3-43          bbmle_1.0.24           numDeriv_2016.8-1.1   
    ## [115] ggbeeswarm_0.6.0
