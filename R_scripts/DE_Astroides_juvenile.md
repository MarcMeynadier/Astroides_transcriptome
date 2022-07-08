DE_Astroides_adult_juveniles
================
Marc Meynadier
6/3/2022

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

    ## Warning: le package 'S4Vectors' a été compilé avec la version R 4.1.3

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

    ## This is vegan 2.6-2

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

    ## Warning: le package 'limma' a été compilé avec la version R 4.1.3

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
```

    ## Note: importing `abundance.h5` is typically faster than `abundance.tsv`

    ## reading in files with read_tsv

    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 
    ## transcripts missing from tx2gene: 1
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
txi2<-tximport(files = files2,type='kallisto',tx2gene = tx2gene)
```

    ## Note: importing `abundance.h5` is typically faster than `abundance.tsv`
    ## reading in files with read_tsv
    ## 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 
    ## transcripts missing from tx2gene: 1
    ## summarizing abundance
    ## summarizing counts
    ## summarizing length

``` r
txi3<-tximport(files = files3,type='kallisto',tx2gene = tx2gene)
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
names(txi2)
```

    ## [1] "abundance"           "counts"              "length"             
    ## [4] "countsFromAbundance"

``` r
names(txi3)
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
    ## TRINITY_DN0_c0_g2                                248.402
    ## TRINITY_DN0_c1_g1                                104.596
    ## TRINITY_DN1_c0_g1                               6103.210
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1704.760
    ##                    abundance_juvenile_rep3_gm_amb_R3AA9_
    ## TRINITY_DN0_c0_g1                              1247.0000
    ## TRINITY_DN0_c0_g2                               140.4570
    ## TRINITY_DN0_c1_g1                                64.3729
    ## TRINITY_DN1_c0_g1                              7743.7300
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1980.4600
    ##                    abundance_juvenile_rep4_gm_amb_R4BA4_
    ## TRINITY_DN0_c0_g1                               1390.000
    ## TRINITY_DN0_c0_g2                                267.593
    ## TRINITY_DN0_c1_g1                                143.146
    ## TRINITY_DN1_c0_g1                               7346.000
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2026.540
    ##                    abundance_juvenile_rep1_gm_low_R1BA8_
    ## TRINITY_DN0_c0_g1                               637.0000
    ## TRINITY_DN0_c0_g2                               137.4210
    ## TRINITY_DN0_c1_g1                                75.0975
    ## TRINITY_DN1_c0_g1                              8165.0000
    ## TRINITY_DN1_c1_g1                                 0.0000
    ## TRINITY_DN10_c0_g1                             1309.2700
    ##                    abundance_juvenile_rep2_gm_low_R2EA3_
    ## TRINITY_DN0_c0_g1                                799.000
    ## TRINITY_DN0_c0_g2                                338.944
    ## TRINITY_DN0_c1_g1                                108.556
    ## TRINITY_DN1_c0_g1                               9692.170
    ## TRINITY_DN1_c1_g1                                  4.000
    ## TRINITY_DN10_c0_g1                              1149.110
    ##                    abundance_juvenile_rep3_gm_low_R3DA1_
    ## TRINITY_DN0_c0_g1                               1099.000
    ## TRINITY_DN0_c0_g2                                332.742
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
    ## TRINITY_DN10_c0_g1                              1575.440
    ##                    abundance_juvenile_rep1_gm_ext_R1EA7_
    ## TRINITY_DN0_c0_g1                                862.000
    ## TRINITY_DN0_c0_g2                                249.717
    ## TRINITY_DN0_c1_g1                                180.206
    ## TRINITY_DN1_c0_g1                               8382.000
    ## TRINITY_DN1_c1_g1                                  2.000
    ## TRINITY_DN10_c0_g1                              1785.440
    ##                    abundance_juvenile_rep2_gm_ext_R2AA2_
    ## TRINITY_DN0_c0_g1                               792.0000
    ## TRINITY_DN0_c0_g2                                98.9913
    ## TRINITY_DN0_c1_g1                               133.4080
    ## TRINITY_DN1_c0_g1                              4871.0000
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1624.2200
    ##                    abundance_juvenile_rep3_gm_ext_R3CA1_
    ## TRINITY_DN0_c0_g1                               1040.000
    ## TRINITY_DN0_c0_g2                                160.000
    ## TRINITY_DN0_c1_g1                                139.177
    ## TRINITY_DN1_c0_g1                               7980.000
    ## TRINITY_DN1_c1_g1                                  2.000
    ## TRINITY_DN10_c0_g1                              1813.330
    ##                    abundance_juvenile_rep4_gm_ext_R4BA15_
    ## TRINITY_DN0_c0_g1                                 707.000
    ## TRINITY_DN0_c0_g2                                 319.875
    ## TRINITY_DN0_c1_g1                                 123.098
    ## TRINITY_DN1_c0_g1                                8087.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               1589.600
    ##                    abundance_juvenile_rep1_sp_amb_R1AA8_
    ## TRINITY_DN0_c0_g1                               1192.000
    ## TRINITY_DN0_c0_g2                                117.034
    ## TRINITY_DN0_c1_g1                                102.447
    ## TRINITY_DN1_c0_g1                               7894.600
    ## TRINITY_DN1_c1_g1                                  6.000
    ## TRINITY_DN10_c0_g1                              1612.620
    ##                    abundance_juvenile_rep2_sp_amb_R2BA7_
    ## TRINITY_DN0_c0_g1                               947.0000
    ## TRINITY_DN0_c0_g2                               224.4740
    ## TRINITY_DN0_c1_g1                                97.6386
    ## TRINITY_DN1_c0_g1                              6051.0000
    ## TRINITY_DN1_c1_g1                                 2.0000
    ## TRINITY_DN10_c0_g1                             1869.8900
    ##                    abundance_juvenile_rep3_sp_amb_R3EA10_
    ## TRINITY_DN0_c0_g1                                1278.000
    ## TRINITY_DN0_c0_g2                                 224.126
    ## TRINITY_DN0_c1_g1                                 120.268
    ## TRINITY_DN1_c0_g1                                6596.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               2252.890
    ##                    abundance_juvenile_rep4_sp_amb_R4AA12_
    ## TRINITY_DN0_c0_g1                               1054.0000
    ## TRINITY_DN0_c0_g2                                168.0370
    ## TRINITY_DN0_c1_g1                                 97.0586
    ## TRINITY_DN1_c0_g1                               7395.0000
    ## TRINITY_DN1_c1_g1                                  0.0000
    ## TRINITY_DN10_c0_g1                              2042.4100
    ##                    abundance_juvenile_rep1_sp_low_R1DA1_
    ## TRINITY_DN0_c0_g1                               923.0000
    ## TRINITY_DN0_c0_g2                               222.0890
    ## TRINITY_DN0_c1_g1                                90.0537
    ## TRINITY_DN1_c0_g1                              6924.0000
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1762.3100
    ##                    abundance_juvenile_rep2_sp_low_R2CA9_
    ## TRINITY_DN0_c0_g1                                924.000
    ## TRINITY_DN0_c0_g2                                208.799
    ## TRINITY_DN0_c1_g1                                132.564
    ## TRINITY_DN1_c0_g1                               6535.180
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1409.990
    ##                    abundance_juvenile_rep3_sp_low_R3BA6_
    ## TRINITY_DN0_c0_g1                               1397.000
    ## TRINITY_DN0_c0_g2                                230.733
    ## TRINITY_DN0_c1_g1                                124.000
    ## TRINITY_DN1_c0_g1                               6544.980
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2067.810
    ##                    abundance_juvenile_rep4_sp_low_R4CA7_
    ## TRINITY_DN0_c0_g1                               1041.000
    ## TRINITY_DN0_c0_g2                                 93.000
    ## TRINITY_DN0_c1_g1                                105.654
    ## TRINITY_DN1_c0_g1                               9480.650
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1832.780

``` r
head(txi2$counts)
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
    ## TRINITY_DN0_c0_g2                                248.402
    ## TRINITY_DN0_c1_g1                                104.596
    ## TRINITY_DN1_c0_g1                               6103.210
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1704.760
    ##                    abundance_juvenile_rep3_gm_amb_R3AA9_
    ## TRINITY_DN0_c0_g1                              1247.0000
    ## TRINITY_DN0_c0_g2                               140.4570
    ## TRINITY_DN0_c1_g1                                64.3729
    ## TRINITY_DN1_c0_g1                              7743.7300
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1980.4600
    ##                    abundance_juvenile_rep4_gm_amb_R4BA4_
    ## TRINITY_DN0_c0_g1                               1390.000
    ## TRINITY_DN0_c0_g2                                267.593
    ## TRINITY_DN0_c1_g1                                143.146
    ## TRINITY_DN1_c0_g1                               7346.000
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2026.540
    ##                    abundance_juvenile_rep1_gm_low_R1BA8_
    ## TRINITY_DN0_c0_g1                               637.0000
    ## TRINITY_DN0_c0_g2                               137.4210
    ## TRINITY_DN0_c1_g1                                75.0975
    ## TRINITY_DN1_c0_g1                              8165.0000
    ## TRINITY_DN1_c1_g1                                 0.0000
    ## TRINITY_DN10_c0_g1                             1309.2700
    ##                    abundance_juvenile_rep2_gm_low_R2EA3_
    ## TRINITY_DN0_c0_g1                                799.000
    ## TRINITY_DN0_c0_g2                                338.944
    ## TRINITY_DN0_c1_g1                                108.556
    ## TRINITY_DN1_c0_g1                               9692.170
    ## TRINITY_DN1_c1_g1                                  4.000
    ## TRINITY_DN10_c0_g1                              1149.110
    ##                    abundance_juvenile_rep3_gm_low_R3DA1_
    ## TRINITY_DN0_c0_g1                               1099.000
    ## TRINITY_DN0_c0_g2                                332.742
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
    ## TRINITY_DN10_c0_g1                              1575.440
    ##                    abundance_juvenile_rep1_sp_amb_R1AA8_
    ## TRINITY_DN0_c0_g1                               1192.000
    ## TRINITY_DN0_c0_g2                                117.034
    ## TRINITY_DN0_c1_g1                                102.447
    ## TRINITY_DN1_c0_g1                               7894.600
    ## TRINITY_DN1_c1_g1                                  6.000
    ## TRINITY_DN10_c0_g1                              1612.620
    ##                    abundance_juvenile_rep2_sp_amb_R2BA7_
    ## TRINITY_DN0_c0_g1                               947.0000
    ## TRINITY_DN0_c0_g2                               224.4740
    ## TRINITY_DN0_c1_g1                                97.6386
    ## TRINITY_DN1_c0_g1                              6051.0000
    ## TRINITY_DN1_c1_g1                                 2.0000
    ## TRINITY_DN10_c0_g1                             1869.8900
    ##                    abundance_juvenile_rep3_sp_amb_R3EA10_
    ## TRINITY_DN0_c0_g1                                1278.000
    ## TRINITY_DN0_c0_g2                                 224.126
    ## TRINITY_DN0_c1_g1                                 120.268
    ## TRINITY_DN1_c0_g1                                6596.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               2252.890
    ##                    abundance_juvenile_rep4_sp_amb_R4AA12_
    ## TRINITY_DN0_c0_g1                               1054.0000
    ## TRINITY_DN0_c0_g2                                168.0370
    ## TRINITY_DN0_c1_g1                                 97.0586
    ## TRINITY_DN1_c0_g1                               7395.0000
    ## TRINITY_DN1_c1_g1                                  0.0000
    ## TRINITY_DN10_c0_g1                              2042.4100
    ##                    abundance_juvenile_rep1_sp_low_R1DA1_
    ## TRINITY_DN0_c0_g1                               923.0000
    ## TRINITY_DN0_c0_g2                               222.0890
    ## TRINITY_DN0_c1_g1                                90.0537
    ## TRINITY_DN1_c0_g1                              6924.0000
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1762.3100
    ##                    abundance_juvenile_rep2_sp_low_R2CA9_
    ## TRINITY_DN0_c0_g1                                924.000
    ## TRINITY_DN0_c0_g2                                208.799
    ## TRINITY_DN0_c1_g1                                132.564
    ## TRINITY_DN1_c0_g1                               6535.180
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1409.990
    ##                    abundance_juvenile_rep3_sp_low_R3BA6_
    ## TRINITY_DN0_c0_g1                               1397.000
    ## TRINITY_DN0_c0_g2                                230.733
    ## TRINITY_DN0_c1_g1                                124.000
    ## TRINITY_DN1_c0_g1                               6544.980
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2067.810
    ##                    abundance_juvenile_rep4_sp_low_R4CA7_
    ## TRINITY_DN0_c0_g1                               1041.000
    ## TRINITY_DN0_c0_g2                                 93.000
    ## TRINITY_DN0_c1_g1                                105.654
    ## TRINITY_DN1_c0_g1                               9480.650
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1832.780

``` r
head(txi3$counts)
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
    ## TRINITY_DN0_c0_g2                                248.402
    ## TRINITY_DN0_c1_g1                                104.596
    ## TRINITY_DN1_c0_g1                               6103.210
    ## TRINITY_DN1_c1_g1                                  0.000
    ## TRINITY_DN10_c0_g1                              1704.760
    ##                    abundance_juvenile_rep3_gm_amb_R3AA9_
    ## TRINITY_DN0_c0_g1                              1247.0000
    ## TRINITY_DN0_c0_g2                               140.4570
    ## TRINITY_DN0_c1_g1                                64.3729
    ## TRINITY_DN1_c0_g1                              7743.7300
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1980.4600
    ##                    abundance_juvenile_rep4_gm_amb_R4BA4_
    ## TRINITY_DN0_c0_g1                               1390.000
    ## TRINITY_DN0_c0_g2                                267.593
    ## TRINITY_DN0_c1_g1                                143.146
    ## TRINITY_DN1_c0_g1                               7346.000
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2026.540
    ##                    abundance_juvenile_rep1_gm_low_R1BA8_
    ## TRINITY_DN0_c0_g1                               637.0000
    ## TRINITY_DN0_c0_g2                               137.4210
    ## TRINITY_DN0_c1_g1                                75.0975
    ## TRINITY_DN1_c0_g1                              8165.0000
    ## TRINITY_DN1_c1_g1                                 0.0000
    ## TRINITY_DN10_c0_g1                             1309.2700
    ##                    abundance_juvenile_rep2_gm_low_R2EA3_
    ## TRINITY_DN0_c0_g1                                799.000
    ## TRINITY_DN0_c0_g2                                338.944
    ## TRINITY_DN0_c1_g1                                108.556
    ## TRINITY_DN1_c0_g1                               9692.170
    ## TRINITY_DN1_c1_g1                                  4.000
    ## TRINITY_DN10_c0_g1                              1149.110
    ##                    abundance_juvenile_rep3_gm_low_R3DA1_
    ## TRINITY_DN0_c0_g1                               1099.000
    ## TRINITY_DN0_c0_g2                                332.742
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
    ## TRINITY_DN10_c0_g1                              1575.440
    ##                    abundance_juvenile_rep1_gm_ext_R1EA7_
    ## TRINITY_DN0_c0_g1                                862.000
    ## TRINITY_DN0_c0_g2                                249.717
    ## TRINITY_DN0_c1_g1                                180.206
    ## TRINITY_DN1_c0_g1                               8382.000
    ## TRINITY_DN1_c1_g1                                  2.000
    ## TRINITY_DN10_c0_g1                              1785.440
    ##                    abundance_juvenile_rep2_gm_ext_R2AA2_
    ## TRINITY_DN0_c0_g1                               792.0000
    ## TRINITY_DN0_c0_g2                                98.9913
    ## TRINITY_DN0_c1_g1                               133.4080
    ## TRINITY_DN1_c0_g1                              4871.0000
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1624.2200
    ##                    abundance_juvenile_rep3_gm_ext_R3CA1_
    ## TRINITY_DN0_c0_g1                               1040.000
    ## TRINITY_DN0_c0_g2                                160.000
    ## TRINITY_DN0_c1_g1                                139.177
    ## TRINITY_DN1_c0_g1                               7980.000
    ## TRINITY_DN1_c1_g1                                  2.000
    ## TRINITY_DN10_c0_g1                              1813.330
    ##                    abundance_juvenile_rep4_gm_ext_R4BA15_
    ## TRINITY_DN0_c0_g1                                 707.000
    ## TRINITY_DN0_c0_g2                                 319.875
    ## TRINITY_DN0_c1_g1                                 123.098
    ## TRINITY_DN1_c0_g1                                8087.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               1589.600
    ##                    abundance_juvenile_rep1_sp_amb_R1AA8_
    ## TRINITY_DN0_c0_g1                               1192.000
    ## TRINITY_DN0_c0_g2                                117.034
    ## TRINITY_DN0_c1_g1                                102.447
    ## TRINITY_DN1_c0_g1                               7894.600
    ## TRINITY_DN1_c1_g1                                  6.000
    ## TRINITY_DN10_c0_g1                              1612.620
    ##                    abundance_juvenile_rep2_sp_amb_R2BA7_
    ## TRINITY_DN0_c0_g1                               947.0000
    ## TRINITY_DN0_c0_g2                               224.4740
    ## TRINITY_DN0_c1_g1                                97.6386
    ## TRINITY_DN1_c0_g1                              6051.0000
    ## TRINITY_DN1_c1_g1                                 2.0000
    ## TRINITY_DN10_c0_g1                             1869.8900
    ##                    abundance_juvenile_rep3_sp_amb_R3EA10_
    ## TRINITY_DN0_c0_g1                                1278.000
    ## TRINITY_DN0_c0_g2                                 224.126
    ## TRINITY_DN0_c1_g1                                 120.268
    ## TRINITY_DN1_c0_g1                                6596.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               2252.890
    ##                    abundance_juvenile_rep4_sp_amb_R4AA12_
    ## TRINITY_DN0_c0_g1                               1054.0000
    ## TRINITY_DN0_c0_g2                                168.0370
    ## TRINITY_DN0_c1_g1                                 97.0586
    ## TRINITY_DN1_c0_g1                               7395.0000
    ## TRINITY_DN1_c1_g1                                  0.0000
    ## TRINITY_DN10_c0_g1                              2042.4100
    ##                    abundance_juvenile_rep1_sp_low_R1DA1_
    ## TRINITY_DN0_c0_g1                               923.0000
    ## TRINITY_DN0_c0_g2                               222.0890
    ## TRINITY_DN0_c1_g1                                90.0537
    ## TRINITY_DN1_c0_g1                              6924.0000
    ## TRINITY_DN1_c1_g1                                 1.0000
    ## TRINITY_DN10_c0_g1                             1762.3100
    ##                    abundance_juvenile_rep2_sp_low_R2CA9_
    ## TRINITY_DN0_c0_g1                                924.000
    ## TRINITY_DN0_c0_g2                                208.799
    ## TRINITY_DN0_c1_g1                                132.564
    ## TRINITY_DN1_c0_g1                               6535.180
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1409.990
    ##                    abundance_juvenile_rep3_sp_low_R3BA6_
    ## TRINITY_DN0_c0_g1                               1397.000
    ## TRINITY_DN0_c0_g2                                230.733
    ## TRINITY_DN0_c1_g1                                124.000
    ## TRINITY_DN1_c0_g1                               6544.980
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              2067.810
    ##                    abundance_juvenile_rep4_sp_low_R4CA7_
    ## TRINITY_DN0_c0_g1                               1041.000
    ## TRINITY_DN0_c0_g2                                 93.000
    ## TRINITY_DN0_c1_g1                                105.654
    ## TRINITY_DN1_c0_g1                               9480.650
    ## TRINITY_DN1_c1_g1                                  1.000
    ## TRINITY_DN10_c0_g1                              1832.780

``` r
head(txiNatSim$counts)
```

    ##                    abundance_juvenile_rep1_gm_low_R1BA8_
    ## TRINITY_DN0_c0_g1                               637.0000
    ## TRINITY_DN0_c0_g2                               137.4210
    ## TRINITY_DN0_c1_g1                                75.0975
    ## TRINITY_DN1_c0_g1                              8165.0000
    ## TRINITY_DN1_c1_g1                                 0.0000
    ## TRINITY_DN10_c0_g1                             1309.2700
    ##                    abundance_juvenile_rep2_gm_low_R2EA3_
    ## TRINITY_DN0_c0_g1                                799.000
    ## TRINITY_DN0_c0_g2                                338.944
    ## TRINITY_DN0_c1_g1                                108.556
    ## TRINITY_DN1_c0_g1                               9692.170
    ## TRINITY_DN1_c1_g1                                  4.000
    ## TRINITY_DN10_c0_g1                              1149.110
    ##                    abundance_juvenile_rep3_gm_low_R3DA1_
    ## TRINITY_DN0_c0_g1                               1099.000
    ## TRINITY_DN0_c0_g2                                332.742
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
    ## TRINITY_DN10_c0_g1                              1575.440
    ##                    abundance_juvenile_rep1_sp_amb_R1AA8_
    ## TRINITY_DN0_c0_g1                               1192.000
    ## TRINITY_DN0_c0_g2                                117.034
    ## TRINITY_DN0_c1_g1                                102.447
    ## TRINITY_DN1_c0_g1                               7894.600
    ## TRINITY_DN1_c1_g1                                  6.000
    ## TRINITY_DN10_c0_g1                              1612.620
    ##                    abundance_juvenile_rep2_sp_amb_R2BA7_
    ## TRINITY_DN0_c0_g1                               947.0000
    ## TRINITY_DN0_c0_g2                               224.4740
    ## TRINITY_DN0_c1_g1                                97.6386
    ## TRINITY_DN1_c0_g1                              6051.0000
    ## TRINITY_DN1_c1_g1                                 2.0000
    ## TRINITY_DN10_c0_g1                             1869.8900
    ##                    abundance_juvenile_rep3_sp_amb_R3EA10_
    ## TRINITY_DN0_c0_g1                                1278.000
    ## TRINITY_DN0_c0_g2                                 224.126
    ## TRINITY_DN0_c1_g1                                 120.268
    ## TRINITY_DN1_c0_g1                                6596.000
    ## TRINITY_DN1_c1_g1                                   2.000
    ## TRINITY_DN10_c0_g1                               2252.890
    ##                    abundance_juvenile_rep4_sp_amb_R4AA12_
    ## TRINITY_DN0_c0_g1                               1054.0000
    ## TRINITY_DN0_c0_g2                                168.0370
    ## TRINITY_DN0_c1_g1                                 97.0586
    ## TRINITY_DN1_c0_g1                               7395.0000
    ## TRINITY_DN1_c1_g1                                  0.0000
    ## TRINITY_DN10_c0_g1                              2042.4100

``` r
dds<-DESeqDataSetFromTximport(txi,colData=samples,design= ~site + pH)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

``` r
dds2<-DESeqDataSetFromTximport(txi2,colData=samples2,design= ~site + pH)
```

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

    ## using counts and average transcript lengths from tximport

``` r
dds3<-DESeqDataSetFromTximport(txi3,colData=samples3,design= ~site_pH)
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
keep <- rowSums(counts(dds2)) >= 10 
dds2 <- dds2[keep,]
keep <- rowSums(counts(dds3)) >= 10 
dds3 <- dds3[keep,]
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
dds2<-DESeq(dds2)
```

    ## estimating size factors
    ## using 'avgTxLength' from assays(dds), correcting for library size
    ## estimating dispersions
    ## gene-wise dispersion estimates
    ## mean-dispersion relationship
    ## final dispersion estimates
    ## fitting model and testing

``` r
dds3<-DESeq(dds3)
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
    ## [2,] "site_SP_vs_GM"            
    ## [3,] "pH_extreme_low_vs_ambient"
    ## [4,] "pH_low_vs_ambient"

``` r
cbind(resultsNames(dds2))
```

    ##      [,1]               
    ## [1,] "Intercept"        
    ## [2,] "site_SP_vs_GM"    
    ## [3,] "pH_low_vs_ambient"

``` r
cbind(resultsNames(dds3))
```

    ##      [,1]                                  
    ## [1,] "Intercept"                           
    ## [2,] "site_pH_GM_extreme_low_vs_GM_ambient"
    ## [3,] "site_pH_GM_low_vs_GM_ambient"        
    ## [4,] "site_pH_SP_ambient_vs_GM_ambient"    
    ## [5,] "site_pH_SP_low_vs_GM_ambient"

``` r
cbind(resultsNames(ddsNatSim))
```

    ##      [,1]                      
    ## [1,] "Intercept"               
    ## [2,] "site_pH_sp_amb_vs_gm_low"

``` r
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
```

    ## 
    ## out of 92052 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 375, 0.41%
    ## LFC < 0 (down)     : 489, 0.53%
    ## outliers [1]       : 564, 0.61%
    ## low counts [2]     : 7110, 7.7%
    ## (mean count < 1)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(sp_VS_gm)
```

    ## 
    ## out of 91025 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 411, 0.45%
    ## LFC < 0 (down)     : 477, 0.52%
    ## outliers [1]       : 732, 0.8%
    ## low counts [2]     : 14009, 15%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(amb_VS_ext)
```

    ## 
    ## out of 92052 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 28, 0.03%
    ## LFC < 0 (down)     : 12, 0.013%
    ## outliers [1]       : 564, 0.61%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(amb_VS_low)
```

    ## 
    ## out of 92052 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 6, 0.0065%
    ## LFC < 0 (down)     : 5, 0.0054%
    ## outliers [1]       : 564, 0.61%
    ## low counts [2]     : 19500, 21%
    ## (mean count < 4)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(low_VS_ext)
```

    ## 
    ## out of 92052 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 37, 0.04%
    ## LFC < 0 (down)     : 10, 0.011%
    ## outliers [1]       : 564, 0.61%
    ## low counts [2]     : 33649, 37%
    ## (mean count < 10)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
summary(sp_amb_VS_gm_low_natSim)
```

    ## 
    ## out of 86188 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 357, 0.41%
    ## LFC < 0 (down)     : 365, 0.42%
    ## outliers [1]       : 663, 0.77%
    ## low counts [2]     : 15027, 17%
    ## (mean count < 4)
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
vsd = vst(dds3,blind=T)

pcaData = plotPCA(vsd, intgroup="site_pH", 
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

pcaData$site_pH = factor(pcaData$site_pH, levels=c("GM_extreme_low","GM_low","GM_ambient","SP_low","SP_ambient"))

ggplot(pcaData, aes(PC1, PC2, colour = site_pH)) + 
  geom_point(size = 5) + theme_bw() + 
  scale_color_manual(values = c("#D55E00","#E69F00","#FFFF00","#0072B2","#56B4E9")) +
  geom_point() +
  theme(text = element_text(size=10),legend.text = element_text(size=10), legend.position = 'bottom',
  axis.text=element_text(size=12), axis.title = element_text(size=12)) +
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

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF","#009E73"),
  stroke_size = 0.4, set_name_size = 4
)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-13.png)<!-- -->

``` r
resOrdered_sp_VS_gm <- sp_VS_gm[order(sp_VS_gm$padj),]
resOrderedDF_sp_VS_gm <- as.data.frame(resOrdered_sp_VS_gm)
resOrderedDF_sp_VS_gm_venn <- filter(resOrderedDF_sp_VS_gm,padj < 0.05)
resOrderedDF_sp_VS_gm_venn <- list(rownames(resOrderedDF_sp_VS_gm_venn))
resOrderedDF_sp_VS_gm_venn <- unlist(resOrderedDF_sp_VS_gm_venn)

x = list('SP VS GM' = resOrderedDF_sp_VS_gm_venn,'ambient VS low' = resOrderedDF_amb_VS_low_venn)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-14.png)<!-- -->

``` r
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

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-15.png)<!-- -->

``` r
x = list('GM ambient VS SP ambient' = resOrderedDF_gm_amb_VS_sp_amb_venn,'GM low VS SP low' = resOrderedDF_gm_low_VS_sp_low_venn)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-16.png)<!-- -->

``` r
x = list('GM ambient VS GM extreme low' = resOrderedDF_gm_amb_VS_gm_ext_venn,'GM low VS GM extreme low' = resOrderedDF_gm_low_VS_gm_ext_venn)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.4, set_name_size = 4
)
```

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-17.png)<!-- -->

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
heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",key.title = "",
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

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-18.png)<!-- -->

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
heatmap.2(assay(vsdCandidate)[topVarGenesVsd,], trace="none",scale="row",keysize=1.15,key.xlab = "",key.title = "",
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

![](DE_Astroides_juvenile_files/figure-gfm/unnamed-chunk-1-19.png)<!-- -->

``` r
# Inferences statistics

count_tab_assay <- assay(vsd)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samples,dist_tab_assay ~ site + pH, method="euclidian")
```

    ## 'adonis' will be deprecated: use 'adonis2' instead

    ## $aov.tab
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
    ## site       1     44036   44036 2.03845 0.10218  0.001 ***
    ## pH         2     41307   20653 0.95604 0.09584  0.702    
    ## Residuals 16    345647   21603         0.80198           
    ## Total     19    430990                 1.00000           
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $call
    ## adonis(formula = dist_tab_assay ~ site + pH, data = samples, 
    ##     method = "euclidian")
    ## 
    ## $coefficients
    ## NULL
    ## 
    ## $coef.sites
    ##                   [,1]      [,2]      [,3]      [,4]       [,5]       [,6]
    ## (Intercept) 198.926574 206.65351 217.17743 215.89223 208.497039 226.552910
    ## site1       -14.875529 -16.19649 -14.29933 -13.88738 -16.503874 -18.791841
    ## pH1         -17.847062 -18.44533 -23.75706 -23.70299  10.722198   9.084261
    ## pH2           9.660025  12.22078  15.95454  13.48571   8.907993  15.085698
    ##                   [,7]       [,8]      [,9]     [,10]      [,11]       [,12]
    ## (Intercept) 209.041641 200.659696 190.45973 190.90378 192.302731 191.6141029
    ## site1       -13.953126 -14.249823  -3.75992  -2.79110  -3.555046  -0.8715483
    ## pH1           3.918304   4.908729  19.74227  19.32897  15.286502  19.7871656
    ## pH2          13.925765   8.758262 -36.95400 -39.08578 -35.416450 -38.5806951
    ##                 [,13]      [,14]      [,15]      [,16]      [,17]     [,18]
    ## (Intercept) 199.40723 191.710004 190.125549 194.501774 193.924116 191.42982
    ## site1        23.59396  18.025028  23.340389  25.043283  16.894388  22.69844
    ## pH1         -12.96658 -11.998993 -11.739353 -11.018634  14.205110  12.24341
    ## pH2          -2.99271  -4.314813  -6.098658  -4.202575  -5.802678  -6.20212
    ##                  [,19]      [,20]
    ## (Intercept) 205.838803 198.552622
    ## site1        20.951255  22.920305
    ## pH1          11.047348  11.531995
    ## pH2          -2.831826  -3.958747
    ## 
    ## $f.perms
    ##              [,1]      [,2]
    ##    [1,] 0.9516879 0.9476042
    ##    [2,] 1.0808890 0.9297500
    ##    [3,] 0.9741367 0.8174119
    ##    [4,] 1.0486950 0.9550258
    ##    [5,] 0.9881476 1.0820896
    ##    [6,] 0.9482389 1.0594358
    ##    [7,] 0.9311127 0.9107699
    ##    [8,] 1.0953274 0.9300092
    ##    [9,] 1.1950588 0.9944616
    ##   [10,] 1.0320729 0.8740032
    ##   [11,] 1.0327612 1.1483867
    ##   [12,] 1.0397246 0.8766347
    ##   [13,] 0.8336192 1.0606975
    ##   [14,] 0.9748327 0.9638110
    ##   [15,] 0.9256616 1.0099543
    ##   [16,] 1.0463858 1.0643185
    ##   [17,] 0.9485353 1.0189120
    ##   [18,] 1.0149341 0.9543095
    ##   [19,] 0.9565701 0.9789075
    ##   [20,] 1.0743337 0.8425766
    ##   [21,] 0.9085044 1.0287970
    ##   [22,] 0.8811003 0.9649398
    ##   [23,] 0.9067693 0.9240314
    ##   [24,] 1.0508342 0.8816370
    ##   [25,] 1.2265210 0.9949182
    ##   [26,] 0.9523699 0.8807300
    ##   [27,] 1.1375865 1.0077521
    ##   [28,] 1.0230226 1.1351635
    ##   [29,] 0.9742621 0.9607683
    ##   [30,] 0.9611679 0.8965303
    ##   [31,] 1.0503035 1.1013395
    ##   [32,] 1.0415804 0.9841225
    ##   [33,] 0.9704101 0.9510623
    ##   [34,] 1.3273079 0.9284692
    ##   [35,] 0.9610453 1.0514351
    ##   [36,] 1.0543670 0.9201261
    ##   [37,] 0.9530997 0.9740945
    ##   [38,] 0.9908791 1.0561604
    ##   [39,] 0.9524826 1.1171648
    ##   [40,] 0.9561825 1.0634289
    ##   [41,] 1.0628180 0.9225348
    ##   [42,] 0.9374950 0.9548553
    ##   [43,] 0.8316432 1.0653747
    ##   [44,] 1.0915431 0.9652400
    ##   [45,] 0.8922082 0.9793117
    ##   [46,] 0.9452508 1.0857560
    ##   [47,] 0.9340551 1.0153388
    ##   [48,] 1.0313938 0.9784427
    ##   [49,] 0.9089878 1.0684537
    ##   [50,] 1.0752856 1.0532112
    ##   [51,] 0.9012482 1.0250665
    ##   [52,] 1.1693852 0.9913434
    ##   [53,] 0.8885602 0.9291534
    ##   [54,] 0.8494440 0.9578945
    ##   [55,] 0.9190479 0.9998622
    ##   [56,] 0.8573901 0.9179396
    ##   [57,] 0.9237808 0.8984190
    ##   [58,] 1.0611718 1.0156524
    ##   [59,] 0.9460687 0.9378688
    ##   [60,] 0.9329668 1.0860199
    ##   [61,] 0.9844313 1.2183064
    ##   [62,] 0.9323044 1.0462717
    ##   [63,] 0.9893975 0.8217382
    ##   [64,] 0.9035456 0.9102639
    ##   [65,] 0.8740750 0.9766547
    ##   [66,] 0.8843452 0.9406785
    ##   [67,] 0.9098584 0.9434014
    ##   [68,] 1.0744702 1.0116768
    ##   [69,] 1.0459711 1.0995359
    ##   [70,] 1.0358074 0.9323464
    ##   [71,] 1.0260079 0.9424996
    ##   [72,] 0.9291413 1.0579175
    ##   [73,] 1.0436885 0.9843946
    ##   [74,] 0.9665994 0.9220790
    ##   [75,] 0.9906915 1.0235595
    ##   [76,] 0.9474759 0.9765887
    ##   [77,] 0.9365410 0.9214430
    ##   [78,] 1.2146854 1.0776422
    ##   [79,] 1.0667771 0.9837678
    ##   [80,] 0.8934621 0.8681782
    ##   [81,] 0.9069163 0.9760872
    ##   [82,] 0.8361716 0.9061049
    ##   [83,] 0.8656589 1.0169959
    ##   [84,] 0.9637113 1.1164849
    ##   [85,] 1.0849369 0.8970877
    ##   [86,] 0.8925880 0.9785844
    ##   [87,] 1.0684152 0.9918147
    ##   [88,] 0.9360186 1.0381621
    ##   [89,] 1.0593831 0.9906397
    ##   [90,] 0.9639209 0.9906754
    ##   [91,] 0.9435574 1.0423353
    ##   [92,] 0.8830364 1.0263311
    ##   [93,] 0.8816130 0.9390963
    ##   [94,] 1.0016354 0.9170113
    ##   [95,] 0.9358280 0.9529471
    ##   [96,] 0.9085794 1.0328698
    ##   [97,] 0.9201454 0.9264430
    ##   [98,] 0.9615679 0.9565998
    ##   [99,] 0.9794183 0.9252620
    ##  [100,] 0.8513961 0.9501409
    ##  [101,] 0.9161522 0.9862594
    ##  [102,] 1.2665004 1.1622714
    ##  [103,] 1.0569152 1.0493712
    ##  [104,] 0.9703009 0.8675086
    ##  [105,] 0.9300115 1.0500949
    ##  [106,] 0.9882356 1.0353780
    ##  [107,] 0.9415193 0.9445089
    ##  [108,] 0.8611300 1.0138464
    ##  [109,] 1.0624616 0.9816578
    ##  [110,] 0.9326690 1.0296138
    ##  [111,] 0.8877981 0.9636344
    ##  [112,] 0.9955497 1.1197861
    ##  [113,] 0.9512694 0.9557355
    ##  [114,] 1.0063424 0.9403490
    ##  [115,] 1.0419477 1.0793239
    ##  [116,] 1.0042191 0.9605746
    ##  [117,] 1.0630702 1.0403437
    ##  [118,] 0.9802592 0.9043380
    ##  [119,] 0.8973256 1.1761107
    ##  [120,] 0.9433187 1.0005895
    ##  [121,] 0.9430708 1.0507860
    ##  [122,] 1.1549133 1.1060948
    ##  [123,] 0.8483400 0.9171884
    ##  [124,] 1.1066580 1.1514425
    ##  [125,] 0.9966794 1.0975313
    ##  [126,] 1.1296866 0.8850158
    ##  [127,] 0.9984191 1.0481172
    ##  [128,] 1.1041060 1.1655372
    ##  [129,] 0.8812161 0.8471810
    ##  [130,] 0.9231683 1.0211792
    ##  [131,] 1.1122108 0.9329310
    ##  [132,] 0.9724693 0.9085358
    ##  [133,] 0.9451572 1.0076142
    ##  [134,] 1.1904256 0.9760745
    ##  [135,] 1.0703756 1.0107560
    ##  [136,] 1.2515976 0.9284481
    ##  [137,] 1.0149548 1.2728851
    ##  [138,] 1.0178340 1.0158337
    ##  [139,] 0.9444794 1.2792897
    ##  [140,] 0.9612025 0.9973263
    ##  [141,] 1.0529372 0.9284820
    ##  [142,] 0.8979892 1.0462851
    ##  [143,] 1.0403030 1.1141563
    ##  [144,] 1.3090512 0.9939130
    ##  [145,] 0.9944517 0.9492487
    ##  [146,] 1.0161572 0.9857168
    ##  [147,] 1.2542699 1.0183158
    ##  [148,] 0.9660703 0.9881718
    ##  [149,] 0.9528026 0.9255066
    ##  [150,] 1.1044657 1.1928199
    ##  [151,] 0.9899997 0.9924941
    ##  [152,] 1.0005405 1.1099386
    ##  [153,] 0.8985879 1.0103259
    ##  [154,] 0.9675337 0.9738760
    ##  [155,] 0.8677240 0.9600711
    ##  [156,] 0.9285444 1.1020269
    ##  [157,] 0.9927750 0.9847290
    ##  [158,] 1.0292344 0.9499826
    ##  [159,] 1.2815490 0.9315164
    ##  [160,] 1.0012943 0.9649089
    ##  [161,] 0.9200913 1.1346906
    ##  [162,] 1.0312854 0.9325040
    ##  [163,] 1.0125677 1.0190190
    ##  [164,] 0.9356137 0.9311247
    ##  [165,] 0.9539257 1.0051419
    ##  [166,] 0.8696432 0.9272905
    ##  [167,] 0.9927107 1.0075426
    ##  [168,] 0.9665290 1.0127094
    ##  [169,] 0.9130548 1.1030374
    ##  [170,] 1.0128282 0.9629710
    ##  [171,] 1.1006641 0.9179222
    ##  [172,] 0.9408827 0.9576439
    ##  [173,] 0.9945543 1.1749412
    ##  [174,] 0.9199377 0.9805231
    ##  [175,] 1.0399692 0.9961794
    ##  [176,] 0.9855966 1.0534868
    ##  [177,] 0.9451787 0.9579659
    ##  [178,] 0.9247618 0.9678264
    ##  [179,] 1.3655911 0.9942555
    ##  [180,] 0.8603250 1.0892125
    ##  [181,] 0.9011878 0.9678908
    ##  [182,] 0.8623911 0.9271080
    ##  [183,] 1.0541613 1.0588933
    ##  [184,] 0.9295792 0.9216912
    ##  [185,] 0.9785360 0.9512958
    ##  [186,] 1.0539525 1.0222948
    ##  [187,] 1.0238153 1.0256221
    ##  [188,] 1.2298380 1.0987009
    ##  [189,] 0.9756305 1.0636508
    ##  [190,] 1.0814643 0.9758931
    ##  [191,] 0.9753048 1.0768599
    ##  [192,] 0.8948500 1.2104797
    ##  [193,] 0.8979144 1.2149613
    ##  [194,] 0.9558216 1.0265597
    ##  [195,] 0.8934243 1.0580379
    ##  [196,] 0.8498966 0.8887639
    ##  [197,] 1.0645683 1.0781303
    ##  [198,] 1.0040381 0.8546556
    ##  [199,] 0.9295563 0.9982873
    ##  [200,] 0.9774181 0.9719211
    ##  [201,] 0.9872557 0.9876453
    ##  [202,] 0.9547601 0.9368210
    ##  [203,] 0.9482949 0.9574216
    ##  [204,] 1.0223791 0.9713785
    ##  [205,] 1.4098436 0.9301248
    ##  [206,] 0.8152974 1.0454319
    ##  [207,] 0.9937383 0.9862034
    ##  [208,] 0.8960185 1.1220389
    ##  [209,] 0.9326885 1.1199483
    ##  [210,] 1.1307667 0.9553283
    ##  [211,] 0.9919399 1.0333198
    ##  [212,] 0.9599149 0.9877369
    ##  [213,] 0.9142394 0.9605914
    ##  [214,] 0.8978045 1.0494344
    ##  [215,] 0.9745460 1.0541427
    ##  [216,] 0.8925890 0.9705865
    ##  [217,] 1.1263363 0.9283372
    ##  [218,] 0.9586402 0.8701187
    ##  [219,] 0.9013502 0.9794762
    ##  [220,] 1.0446162 0.8256038
    ##  [221,] 0.9698642 1.3021525
    ##  [222,] 0.9808019 0.9695345
    ##  [223,] 0.9387631 1.0794221
    ##  [224,] 0.9447365 1.0824458
    ##  [225,] 1.0169738 0.9567325
    ##  [226,] 0.9755390 1.0883128
    ##  [227,] 0.8814784 0.9556411
    ##  [228,] 0.8887595 1.0384977
    ##  [229,] 1.0861580 1.0562795
    ##  [230,] 1.0089766 1.1135344
    ##  [231,] 0.9624652 1.1796565
    ##  [232,] 0.9655168 0.9546464
    ##  [233,] 1.0628716 1.0197679
    ##  [234,] 0.9032153 0.9025944
    ##  [235,] 0.9168802 1.0745580
    ##  [236,] 0.9043512 0.8719942
    ##  [237,] 1.0390570 1.1033363
    ##  [238,] 1.0845082 1.0149618
    ##  [239,] 1.0298671 1.1295407
    ##  [240,] 0.9933461 1.0078225
    ##  [241,] 0.9748950 0.9462064
    ##  [242,] 1.0169022 1.0717384
    ##  [243,] 1.0204822 0.8772450
    ##  [244,] 0.7585008 0.9518388
    ##  [245,] 0.9891426 0.9897004
    ##  [246,] 0.8992989 1.0934382
    ##  [247,] 0.9093751 0.9047791
    ##  [248,] 0.8142181 0.9464306
    ##  [249,] 0.9380409 1.0284864
    ##  [250,] 1.0647402 0.9993276
    ##  [251,] 1.1168883 1.1935505
    ##  [252,] 1.2280055 1.3305686
    ##  [253,] 1.2490182 1.0311575
    ##  [254,] 1.0959870 1.0378575
    ##  [255,] 0.9915655 0.9568443
    ##  [256,] 0.9158228 1.0037448
    ##  [257,] 0.9734666 1.0048238
    ##  [258,] 0.9181780 1.1843042
    ##  [259,] 0.9271836 0.9933997
    ##  [260,] 0.9992013 1.0393395
    ##  [261,] 0.9001638 1.0200610
    ##  [262,] 1.2233757 1.0583797
    ##  [263,] 1.2174817 1.0197185
    ##  [264,] 0.9714189 1.1450110
    ##  [265,] 0.9432957 1.0254052
    ##  [266,] 1.1011683 1.2568885
    ##  [267,] 0.9619400 1.0762390
    ##  [268,] 0.9213397 0.8164879
    ##  [269,] 0.9083344 1.1176096
    ##  [270,] 0.9837952 0.8844876
    ##  [271,] 0.9026255 0.9297246
    ##  [272,] 1.1994090 0.9909320
    ##  [273,] 0.8880732 0.9084289
    ##  [274,] 1.1737008 1.0053070
    ##  [275,] 0.9925042 1.0064478
    ##  [276,] 1.0923516 1.1209938
    ##  [277,] 1.0590908 1.1727430
    ##  [278,] 0.9969909 0.9853524
    ##  [279,] 0.8986293 1.0744634
    ##  [280,] 0.9694456 0.9133342
    ##  [281,] 0.8903653 0.9677305
    ##  [282,] 0.9352970 0.9322529
    ##  [283,] 1.0058451 0.9634309
    ##  [284,] 0.9191757 0.9219776
    ##  [285,] 1.0514468 1.0268198
    ##  [286,] 1.1466657 1.0610005
    ##  [287,] 0.9481831 1.1020555
    ##  [288,] 1.1545064 0.9045957
    ##  [289,] 0.8969786 1.0167864
    ##  [290,] 0.9187032 1.0065992
    ##  [291,] 1.0257439 0.9173795
    ##  [292,] 0.8695014 0.8631387
    ##  [293,] 1.0329187 0.9481032
    ##  [294,] 0.9310022 0.9577448
    ##  [295,] 0.9321043 0.9288263
    ##  [296,] 1.1394431 1.0265338
    ##  [297,] 1.2670931 1.0323069
    ##  [298,] 0.9541813 1.0096568
    ##  [299,] 0.8951266 0.9524927
    ##  [300,] 0.9824086 0.9367824
    ##  [301,] 0.9592033 0.9676976
    ##  [302,] 1.0720503 0.9055068
    ##  [303,] 0.9774958 1.0904567
    ##  [304,] 1.1084125 0.8709728
    ##  [305,] 0.9816820 0.9561637
    ##  [306,] 1.1297377 1.0524301
    ##  [307,] 0.9826178 1.0466604
    ##  [308,] 0.9894126 0.9005539
    ##  [309,] 0.8853000 1.0142892
    ##  [310,] 1.0652864 1.1085359
    ##  [311,] 1.0160197 0.9175538
    ##  [312,] 0.9648593 0.9579289
    ##  [313,] 0.8915406 0.9475392
    ##  [314,] 1.1663283 1.3899525
    ##  [315,] 0.8167429 1.0021658
    ##  [316,] 0.8739392 1.0501630
    ##  [317,] 0.9628436 0.9150449
    ##  [318,] 0.8949373 0.9655393
    ##  [319,] 0.9022887 1.1081534
    ##  [320,] 0.9696477 1.0798995
    ##  [321,] 1.0249723 0.9001112
    ##  [322,] 1.0177916 0.9836460
    ##  [323,] 0.9305365 1.1451270
    ##  [324,] 0.9436200 0.9568972
    ##  [325,] 0.9646388 0.9636723
    ##  [326,] 0.9882983 1.0478481
    ##  [327,] 1.0030678 1.0262025
    ##  [328,] 1.3511250 1.0034015
    ##  [329,] 0.9146398 0.9858739
    ##  [330,] 1.0436601 0.9815151
    ##  [331,] 1.0294538 0.9179150
    ##  [332,] 1.0992477 1.0090088
    ##  [333,] 1.0637982 1.0268877
    ##  [334,] 0.9586225 0.9876977
    ##  [335,] 0.9335908 1.0277977
    ##  [336,] 1.1930027 1.0467627
    ##  [337,] 0.8863882 1.0087918
    ##  [338,] 0.9429276 1.0525280
    ##  [339,] 1.3008348 0.9176116
    ##  [340,] 0.9794887 0.9895768
    ##  [341,] 1.1161836 1.0088637
    ##  [342,] 1.0564698 1.0889235
    ##  [343,] 1.0102053 0.9929559
    ##  [344,] 0.9143436 0.9759290
    ##  [345,] 0.9182971 0.9487229
    ##  [346,] 0.9334457 1.0063418
    ##  [347,] 1.1769585 1.0563255
    ##  [348,] 1.0724010 0.8384787
    ##  [349,] 1.0865048 1.0047864
    ##  [350,] 1.1210381 0.9351483
    ##  [351,] 0.8815754 1.1392637
    ##  [352,] 1.2082911 0.9893539
    ##  [353,] 1.0494189 1.0215476
    ##  [354,] 1.1565331 1.0697436
    ##  [355,] 1.1248761 1.1244551
    ##  [356,] 1.0851786 1.1697258
    ##  [357,] 1.0119783 0.8915924
    ##  [358,] 0.9586190 1.0401229
    ##  [359,] 0.9309415 1.0397890
    ##  [360,] 1.0033010 0.9226221
    ##  [361,] 1.0346201 1.0104869
    ##  [362,] 0.9446074 0.8469198
    ##  [363,] 0.9717266 1.0253171
    ##  [364,] 1.0843850 1.2678506
    ##  [365,] 1.2810957 0.9320389
    ##  [366,] 1.0720887 1.0551225
    ##  [367,] 1.0036774 0.9667798
    ##  [368,] 1.0317880 0.9939451
    ##  [369,] 0.9940643 1.0298544
    ##  [370,] 0.8208764 0.8905931
    ##  [371,] 1.0775896 0.8752881
    ##  [372,] 0.9977008 1.1410888
    ##  [373,] 1.0638723 0.9722452
    ##  [374,] 1.0791176 0.9803096
    ##  [375,] 0.8721724 0.9311879
    ##  [376,] 0.9403116 0.9598412
    ##  [377,] 1.0037212 0.9931699
    ##  [378,] 0.9614814 1.1476724
    ##  [379,] 1.0036470 0.9246174
    ##  [380,] 0.9211553 0.8966906
    ##  [381,] 0.9687388 1.1183771
    ##  [382,] 1.0162389 1.0658169
    ##  [383,] 1.0189531 0.9833627
    ##  [384,] 0.9417025 0.9486368
    ##  [385,] 0.9815062 1.0005592
    ##  [386,] 1.0160574 1.0766129
    ##  [387,] 0.9433705 0.9776097
    ##  [388,] 0.9653604 0.9793344
    ##  [389,] 1.1920240 1.2515163
    ##  [390,] 0.8304443 1.0135756
    ##  [391,] 0.9304933 0.9863729
    ##  [392,] 1.0553089 0.9857470
    ##  [393,] 0.8760794 1.0398102
    ##  [394,] 1.0649574 0.9845138
    ##  [395,] 0.9834782 0.8616849
    ##  [396,] 1.0890921 0.9491014
    ##  [397,] 1.4171030 0.9090884
    ##  [398,] 1.1525807 0.9739665
    ##  [399,] 0.8850236 0.9032334
    ##  [400,] 1.0323910 0.9614016
    ##  [401,] 1.0965715 0.9626722
    ##  [402,] 0.8701504 0.9347369
    ##  [403,] 1.2673534 1.0150830
    ##  [404,] 1.0660134 0.9719650
    ##  [405,] 0.8907304 1.0623102
    ##  [406,] 0.9118866 0.9870981
    ##  [407,] 1.0413262 0.9102212
    ##  [408,] 1.0261389 1.1207252
    ##  [409,] 0.9781316 0.9577756
    ##  [410,] 1.0297550 1.0051338
    ##  [411,] 0.8677598 0.9756820
    ##  [412,] 0.9188428 1.2114537
    ##  [413,] 0.9245420 0.9557546
    ##  [414,] 0.9736537 0.8463995
    ##  [415,] 1.2658062 0.9760414
    ##  [416,] 0.9900629 0.9134774
    ##  [417,] 1.1802232 0.9575129
    ##  [418,] 1.4472059 1.0998165
    ##  [419,] 1.0527530 1.0313816
    ##  [420,] 1.0139173 0.9293187
    ##  [421,] 0.9291254 0.9813350
    ##  [422,] 0.9484469 1.0890701
    ##  [423,] 0.8843875 0.8983282
    ##  [424,] 0.9585336 1.1322630
    ##  [425,] 1.0865589 1.1592371
    ##  [426,] 0.9787646 0.9233857
    ##  [427,] 0.8596347 0.9333618
    ##  [428,] 0.9049562 1.2543657
    ##  [429,] 0.8533966 1.1794069
    ##  [430,] 0.9810478 1.0986578
    ##  [431,] 1.0878432 1.0668184
    ##  [432,] 1.0974502 0.9856389
    ##  [433,] 0.8879490 1.0357960
    ##  [434,] 1.0745295 0.8677786
    ##  [435,] 1.1706665 1.0539166
    ##  [436,] 0.9281805 1.0692103
    ##  [437,] 0.8474770 0.8695430
    ##  [438,] 0.9294915 0.9874588
    ##  [439,] 0.8874206 1.0107822
    ##  [440,] 1.0391712 0.9158037
    ##  [441,] 1.1117269 0.8833874
    ##  [442,] 0.9035719 1.0596568
    ##  [443,] 1.0702118 1.0454173
    ##  [444,] 0.9692019 1.0142375
    ##  [445,] 0.8732613 1.1162306
    ##  [446,] 0.9169230 0.9271558
    ##  [447,] 0.9490631 1.0083467
    ##  [448,] 1.0099678 1.2181813
    ##  [449,] 0.9658571 0.9606407
    ##  [450,] 0.9050009 0.9715794
    ##  [451,] 0.9487961 0.9882657
    ##  [452,] 1.0658768 1.0538910
    ##  [453,] 1.0611754 0.9305258
    ##  [454,] 0.9433150 0.8781552
    ##  [455,] 1.2297942 0.9536637
    ##  [456,] 0.9477040 0.9623090
    ##  [457,] 1.0146053 0.9019536
    ##  [458,] 1.2718939 1.0640078
    ##  [459,] 0.8631994 0.9301857
    ##  [460,] 1.0624423 1.0588427
    ##  [461,] 1.0898606 0.8808056
    ##  [462,] 1.1515533 1.0410671
    ##  [463,] 0.9325980 1.0053925
    ##  [464,] 0.9993796 0.9870008
    ##  [465,] 1.0492332 0.9302327
    ##  [466,] 1.7032851 0.8835767
    ##  [467,] 0.9682547 0.9798052
    ##  [468,] 1.0043800 0.8989374
    ##  [469,] 1.0600765 1.1701096
    ##  [470,] 0.9278680 0.9676641
    ##  [471,] 1.0095784 0.9727566
    ##  [472,] 0.9662402 1.0457577
    ##  [473,] 1.0118052 1.0587064
    ##  [474,] 0.9732758 0.9857306
    ##  [475,] 1.1080321 0.9028342
    ##  [476,] 1.0701221 0.9794324
    ##  [477,] 1.1344230 0.9648518
    ##  [478,] 1.0682365 1.1783096
    ##  [479,] 0.9878020 1.0809512
    ##  [480,] 0.9221561 0.9658237
    ##  [481,] 0.9510111 0.9916853
    ##  [482,] 1.0678319 0.9498057
    ##  [483,] 0.9380473 1.0815186
    ##  [484,] 1.0453056 0.9852111
    ##  [485,] 1.0738121 0.9408871
    ##  [486,] 0.8890212 0.9095637
    ##  [487,] 0.9833364 1.0011836
    ##  [488,] 0.8930341 1.0402434
    ##  [489,] 0.9537870 0.9331263
    ##  [490,] 1.0557173 1.0684198
    ##  [491,] 0.9101218 1.0618186
    ##  [492,] 1.1002675 1.0020985
    ##  [493,] 1.1430443 0.9098083
    ##  [494,] 0.8848348 0.9562420
    ##  [495,] 1.0379310 0.9865085
    ##  [496,] 0.9295714 1.0343884
    ##  [497,] 1.2141636 0.9587276
    ##  [498,] 0.9179671 1.1487301
    ##  [499,] 0.9646884 1.1113060
    ##  [500,] 1.1255429 1.0771424
    ##  [501,] 0.9793116 0.9917888
    ##  [502,] 1.0699131 1.0058717
    ##  [503,] 0.9363747 1.1536604
    ##  [504,] 0.8697873 0.9719932
    ##  [505,] 0.9583134 0.9673640
    ##  [506,] 0.9297070 0.9126165
    ##  [507,] 0.9859904 1.1296245
    ##  [508,] 1.1138569 0.9862765
    ##  [509,] 0.8893814 0.9706607
    ##  [510,] 1.1009460 0.9949776
    ##  [511,] 1.1831634 0.9452306
    ##  [512,] 0.9875541 1.0707423
    ##  [513,] 0.9399232 1.1005606
    ##  [514,] 1.0925456 0.8956515
    ##  [515,] 0.8570631 0.9429245
    ##  [516,] 0.8593562 0.9907468
    ##  [517,] 0.8385179 0.8990915
    ##  [518,] 0.7789905 0.9205337
    ##  [519,] 1.0899053 1.0841414
    ##  [520,] 0.9137819 0.8975929
    ##  [521,] 0.8993032 0.9722361
    ##  [522,] 1.0679333 0.8993065
    ##  [523,] 0.9916777 1.1762215
    ##  [524,] 0.9363660 1.0263505
    ##  [525,] 0.8778006 0.9198125
    ##  [526,] 0.9447662 0.9510273
    ##  [527,] 1.0578857 1.0814430
    ##  [528,] 1.0495873 0.9753104
    ##  [529,] 0.8536878 0.9665043
    ##  [530,] 0.9608918 1.0483532
    ##  [531,] 0.9422988 0.9590613
    ##  [532,] 0.9620598 0.9648774
    ##  [533,] 1.0339008 1.2511179
    ##  [534,] 1.0089427 0.9216347
    ##  [535,] 0.9257235 0.9359025
    ##  [536,] 0.8925437 0.9385925
    ##  [537,] 1.4804030 1.1793788
    ##  [538,] 0.8812231 0.8947876
    ##  [539,] 0.9261113 0.9734655
    ##  [540,] 0.8314754 0.9402905
    ##  [541,] 0.8867225 1.0804105
    ##  [542,] 0.9594908 0.9651106
    ##  [543,] 0.9893894 0.9013325
    ##  [544,] 0.9431243 1.0118301
    ##  [545,] 0.9025912 0.9250213
    ##  [546,] 0.9490997 0.9932745
    ##  [547,] 1.1817316 1.0406869
    ##  [548,] 0.9496019 0.9729519
    ##  [549,] 1.0476731 1.0627245
    ##  [550,] 1.1350666 0.9427597
    ##  [551,] 0.9940395 1.0583664
    ##  [552,] 1.0871969 0.9191430
    ##  [553,] 1.0092592 1.1962654
    ##  [554,] 1.0804206 1.0350850
    ##  [555,] 0.9981170 1.2282377
    ##  [556,] 0.9587292 0.9008749
    ##  [557,] 0.9068592 0.9422420
    ##  [558,] 0.8544436 0.9710882
    ##  [559,] 0.9935661 0.9874261
    ##  [560,] 0.8881690 0.9387619
    ##  [561,] 0.8889203 1.2217934
    ##  [562,] 1.3123162 1.0380727
    ##  [563,] 0.9619465 0.9871961
    ##  [564,] 0.8830480 1.0824226
    ##  [565,] 0.9883447 0.8797818
    ##  [566,] 0.9202096 1.0519910
    ##  [567,] 0.9820140 0.9570321
    ##  [568,] 0.9495070 1.0104566
    ##  [569,] 0.8324041 1.0151559
    ##  [570,] 0.9312839 0.9734153
    ##  [571,] 1.0924537 1.3317102
    ##  [572,] 0.9275596 1.0816086
    ##  [573,] 0.9587392 0.9120649
    ##  [574,] 0.9118020 1.0290326
    ##  [575,] 1.1129144 1.1322665
    ##  [576,] 1.1513558 0.9599695
    ##  [577,] 1.0138529 1.0853719
    ##  [578,] 0.9089388 0.9630838
    ##  [579,] 1.0894944 0.9874770
    ##  [580,] 1.0520194 0.9852253
    ##  [581,] 0.9466582 0.8645778
    ##  [582,] 1.0520750 0.9801808
    ##  [583,] 0.8918495 0.9715320
    ##  [584,] 0.8947624 0.9245415
    ##  [585,] 1.0552395 1.1618048
    ##  [586,] 1.1061513 1.0234481
    ##  [587,] 1.0643879 1.0436445
    ##  [588,] 1.1494338 0.9512730
    ##  [589,] 1.0013583 0.9895670
    ##  [590,] 1.0600389 1.0391729
    ##  [591,] 1.2316125 0.8851027
    ##  [592,] 1.1210163 0.9437511
    ##  [593,] 1.0563751 1.0858574
    ##  [594,] 0.9686662 1.0483167
    ##  [595,] 0.9861709 0.8980835
    ##  [596,] 0.9708260 1.0536757
    ##  [597,] 0.9838448 0.9422660
    ##  [598,] 1.0556416 0.9269438
    ##  [599,] 0.8611618 0.9648074
    ##  [600,] 0.9801707 1.0837923
    ##  [601,] 0.9168663 0.8694195
    ##  [602,] 1.0070694 1.0234632
    ##  [603,] 1.0166371 1.1211726
    ##  [604,] 0.8251929 1.3565450
    ##  [605,] 0.9449378 1.1382895
    ##  [606,] 0.9277944 1.1365085
    ##  [607,] 1.0117428 0.9452974
    ##  [608,] 0.8802076 1.1046263
    ##  [609,] 1.1034648 0.9009599
    ##  [610,] 0.8959453 1.0376608
    ##  [611,] 1.0489865 1.1028083
    ##  [612,] 0.9806805 1.1525639
    ##  [613,] 0.8967985 0.9339937
    ##  [614,] 0.9697821 0.9910845
    ##  [615,] 1.1247083 1.0506599
    ##  [616,] 0.9357230 0.9538085
    ##  [617,] 1.0381428 0.9184952
    ##  [618,] 1.3165590 1.0877634
    ##  [619,] 0.8837473 0.9413344
    ##  [620,] 0.8594025 0.9747182
    ##  [621,] 0.8789531 1.0175754
    ##  [622,] 0.9210388 0.9218199
    ##  [623,] 1.0520058 0.9156256
    ##  [624,] 0.9758825 0.8831627
    ##  [625,] 0.9226063 1.0808735
    ##  [626,] 1.0415769 0.9735231
    ##  [627,] 1.0373638 1.0720615
    ##  [628,] 0.9017616 1.1250305
    ##  [629,] 0.8965052 0.9336388
    ##  [630,] 0.9413270 1.0176230
    ##  [631,] 0.9332397 0.9398306
    ##  [632,] 1.0948691 1.1090143
    ##  [633,] 0.8499698 0.9377142
    ##  [634,] 0.9993787 1.1417166
    ##  [635,] 1.2597128 1.1657928
    ##  [636,] 0.9039671 0.9467803
    ##  [637,] 1.1258319 1.0212633
    ##  [638,] 0.8736233 0.9513735
    ##  [639,] 0.9862951 1.0133003
    ##  [640,] 0.9775145 1.0471786
    ##  [641,] 1.0497445 0.9351528
    ##  [642,] 0.9374310 1.0024488
    ##  [643,] 0.9904957 1.0079682
    ##  [644,] 0.9244591 1.0745232
    ##  [645,] 0.9995673 0.9985839
    ##  [646,] 1.0390829 0.9842959
    ##  [647,] 1.0041791 0.9278270
    ##  [648,] 1.0451179 1.1517550
    ##  [649,] 0.8271263 1.0596002
    ##  [650,] 0.9877605 0.9289780
    ##  [651,] 1.0757274 0.9417481
    ##  [652,] 0.8455921 1.0226875
    ##  [653,] 0.8660341 0.9792440
    ##  [654,] 0.9868165 1.0309203
    ##  [655,] 0.9005500 0.9604649
    ##  [656,] 0.9736478 1.0505283
    ##  [657,] 0.9859345 0.9912290
    ##  [658,] 1.0859051 1.0532546
    ##  [659,] 0.8200611 0.9794212
    ##  [660,] 0.8802771 1.0659035
    ##  [661,] 0.9086290 0.9220145
    ##  [662,] 0.9981111 1.1809053
    ##  [663,] 1.0048865 0.9661986
    ##  [664,] 0.9656606 0.9793838
    ##  [665,] 1.0578412 0.9999421
    ##  [666,] 0.8804063 0.8888247
    ##  [667,] 1.0152991 0.9738190
    ##  [668,] 1.0007525 1.3413440
    ##  [669,] 0.9694407 1.0902566
    ##  [670,] 0.9234342 1.0325576
    ##  [671,] 0.9523020 1.0355558
    ##  [672,] 0.9693385 1.0331262
    ##  [673,] 1.0970795 0.9346151
    ##  [674,] 1.0409686 0.9412442
    ##  [675,] 0.8945809 0.9380701
    ##  [676,] 0.9199525 0.8745554
    ##  [677,] 0.8503571 1.0405580
    ##  [678,] 1.0096439 0.9551017
    ##  [679,] 1.1562508 0.9900130
    ##  [680,] 1.1617323 1.0727065
    ##  [681,] 0.9732124 0.9265899
    ##  [682,] 1.0761540 1.0999200
    ##  [683,] 0.9278914 0.9175170
    ##  [684,] 0.9117703 0.9995125
    ##  [685,] 0.9805033 1.0414167
    ##  [686,] 1.2126084 0.9325636
    ##  [687,] 0.9501842 1.2233099
    ##  [688,] 0.9722677 1.0362696
    ##  [689,] 1.0494575 0.8851249
    ##  [690,] 0.9618906 0.9699513
    ##  [691,] 1.1856026 1.0258188
    ##  [692,] 1.0413368 1.0621685
    ##  [693,] 0.8894288 1.2451112
    ##  [694,] 0.9115642 0.9338754
    ##  [695,] 0.9091293 1.0027523
    ##  [696,] 0.8861404 1.0100237
    ##  [697,] 0.8922589 0.9927463
    ##  [698,] 1.0004254 0.9055483
    ##  [699,] 0.9920638 1.1052057
    ##  [700,] 0.9080584 1.1540695
    ##  [701,] 0.9386903 1.2292781
    ##  [702,] 0.9617605 0.9409547
    ##  [703,] 0.9150979 0.8938339
    ##  [704,] 1.1429993 1.1085209
    ##  [705,] 1.0104822 0.9849159
    ##  [706,] 0.8862418 1.0370919
    ##  [707,] 1.0342507 0.9706316
    ##  [708,] 1.2911092 1.1199104
    ##  [709,] 1.3363293 1.0917014
    ##  [710,] 1.0704444 0.9962639
    ##  [711,] 0.9658100 0.9346557
    ##  [712,] 1.1087349 0.9771996
    ##  [713,] 1.0866758 1.0503270
    ##  [714,] 0.9230082 1.2190939
    ##  [715,] 0.8845154 0.9700848
    ##  [716,] 1.0564496 1.0943969
    ##  [717,] 0.9559669 1.0033834
    ##  [718,] 1.0181213 0.9366246
    ##  [719,] 0.9214788 1.0318808
    ##  [720,] 0.8821588 0.9688265
    ##  [721,] 1.0077497 0.9515225
    ##  [722,] 0.9485227 0.9931619
    ##  [723,] 0.8910843 1.0109837
    ##  [724,] 0.9696691 0.8685711
    ##  [725,] 0.8525484 0.9681534
    ##  [726,] 0.9229862 1.1865227
    ##  [727,] 1.2917609 1.1355859
    ##  [728,] 1.3874412 0.9335681
    ##  [729,] 0.8969782 0.8992789
    ##  [730,] 1.3497445 1.0224886
    ##  [731,] 0.9181073 1.2060087
    ##  [732,] 1.0589590 0.9058697
    ##  [733,] 1.1074753 1.1646254
    ##  [734,] 1.0539503 1.0901576
    ##  [735,] 0.8970861 1.1039548
    ##  [736,] 1.1915814 1.1147384
    ##  [737,] 1.0066776 0.9309357
    ##  [738,] 0.9147753 0.9188979
    ##  [739,] 1.0804668 0.9637142
    ##  [740,] 0.9039391 0.9848244
    ##  [741,] 0.9910254 0.9935545
    ##  [742,] 0.9747513 0.9866028
    ##  [743,] 0.9988568 0.9812080
    ##  [744,] 1.0579381 1.0380666
    ##  [745,] 0.9251569 1.0662542
    ##  [746,] 0.8533492 0.9465242
    ##  [747,] 1.0877673 0.9665338
    ##  [748,] 0.8997925 0.9701978
    ##  [749,] 0.9197366 0.9084637
    ##  [750,] 0.8870655 0.9814344
    ##  [751,] 1.1888136 1.0330570
    ##  [752,] 0.9437812 1.0434390
    ##  [753,] 0.9322896 1.0660464
    ##  [754,] 0.9579841 0.8981544
    ##  [755,] 0.8749825 0.9560918
    ##  [756,] 0.9891380 0.9513605
    ##  [757,] 0.9166495 1.0213358
    ##  [758,] 0.9761058 1.1602600
    ##  [759,] 1.0024942 0.9913717
    ##  [760,] 1.0502105 0.9428310
    ##  [761,] 0.9888023 0.8881513
    ##  [762,] 1.0116480 0.9670997
    ##  [763,] 0.9404552 1.0254075
    ##  [764,] 1.1200230 0.9817342
    ##  [765,] 0.9569715 0.9660430
    ##  [766,] 0.9919013 0.9181669
    ##  [767,] 0.9924315 0.9657745
    ##  [768,] 1.0940953 1.0099859
    ##  [769,] 0.8695614 0.9806997
    ##  [770,] 1.0870836 0.9577779
    ##  [771,] 1.0293376 1.1120534
    ##  [772,] 1.0525855 0.9991752
    ##  [773,] 0.9925630 0.9455654
    ##  [774,] 0.8604149 1.0058084
    ##  [775,] 0.9359922 1.1102879
    ##  [776,] 1.0387449 1.1831936
    ##  [777,] 0.9610593 0.9275610
    ##  [778,] 0.9849316 0.9793184
    ##  [779,] 0.8876562 1.0031806
    ##  [780,] 0.9603695 1.0171780
    ##  [781,] 0.9200490 0.9305872
    ##  [782,] 1.0091345 0.9197875
    ##  [783,] 0.9600820 0.9826453
    ##  [784,] 1.3470137 1.0010843
    ##  [785,] 1.5015181 0.9164087
    ##  [786,] 0.9664301 1.0794031
    ##  [787,] 1.0066863 0.9152749
    ##  [788,] 0.9824586 1.0317481
    ##  [789,] 0.9083285 0.9357758
    ##  [790,] 0.9746473 0.9958215
    ##  [791,] 1.0388730 0.9734677
    ##  [792,] 1.0179926 0.9649644
    ##  [793,] 0.9432301 0.9903631
    ##  [794,] 0.8709809 0.9108745
    ##  [795,] 1.1070339 0.9318307
    ##  [796,] 0.8849536 0.9633405
    ##  [797,] 1.2896332 0.9923264
    ##  [798,] 1.2305195 1.0177555
    ##  [799,] 0.9740470 1.1513562
    ##  [800,] 0.9412225 1.0606522
    ##  [801,] 1.2487511 0.9582119
    ##  [802,] 1.2235964 1.0261942
    ##  [803,] 0.9956868 0.9306069
    ##  [804,] 0.9271153 1.0741684
    ##  [805,] 0.8577422 1.0212767
    ##  [806,] 0.8558893 0.9026296
    ##  [807,] 1.0068883 1.1551806
    ##  [808,] 0.8549076 0.9711240
    ##  [809,] 1.0953275 1.0331087
    ##  [810,] 1.0447325 0.9607291
    ##  [811,] 0.9447774 1.0985513
    ##  [812,] 0.9237597 1.1228007
    ##  [813,] 1.0585920 1.0672669
    ##  [814,] 0.9375708 1.0145470
    ##  [815,] 1.1679714 0.9447646
    ##  [816,] 0.9161756 0.8736356
    ##  [817,] 1.0192248 0.9747502
    ##  [818,] 0.9849362 1.0740468
    ##  [819,] 1.0685360 0.9553764
    ##  [820,] 0.9142901 0.9750397
    ##  [821,] 1.0118575 0.9687207
    ##  [822,] 1.2582572 1.0472255
    ##  [823,] 0.9049690 0.9441949
    ##  [824,] 1.1094440 1.0745066
    ##  [825,] 0.9479764 0.8727970
    ##  [826,] 1.1234901 0.9414275
    ##  [827,] 1.0238495 0.8915804
    ##  [828,] 1.0332148 1.0572749
    ##  [829,] 0.8723742 0.9727752
    ##  [830,] 1.0075436 1.0906831
    ##  [831,] 0.8702748 1.1032297
    ##  [832,] 0.8265306 0.9251848
    ##  [833,] 0.8446518 0.9305582
    ##  [834,] 1.1318528 1.0004667
    ##  [835,] 1.0030404 0.9827900
    ##  [836,] 1.0389845 1.1385108
    ##  [837,] 0.9149757 1.0388158
    ##  [838,] 0.9455616 1.0293824
    ##  [839,] 1.0582248 1.1261815
    ##  [840,] 1.0175090 0.9928332
    ##  [841,] 0.8334401 0.9266476
    ##  [842,] 1.0837420 1.0867170
    ##  [843,] 0.9840596 0.9598803
    ##  [844,] 1.1409877 0.9960675
    ##  [845,] 0.9031451 0.9419264
    ##  [846,] 0.8975947 0.9323864
    ##  [847,] 0.9374763 0.8599742
    ##  [848,] 0.9243672 1.0304791
    ##  [849,] 0.9790084 0.9458036
    ##  [850,] 1.0621827 0.9555611
    ##  [851,] 1.0041214 0.9832659
    ##  [852,] 1.1696770 0.9326988
    ##  [853,] 1.2419513 1.0505524
    ##  [854,] 0.9801525 1.0739073
    ##  [855,] 0.8391896 0.9554668
    ##  [856,] 0.9370315 1.0025725
    ##  [857,] 0.9078069 0.9384729
    ##  [858,] 0.8844610 0.9867081
    ##  [859,] 0.9843606 0.9017682
    ##  [860,] 1.1154681 0.9302740
    ##  [861,] 1.0106415 1.0727220
    ##  [862,] 1.0318832 1.1702362
    ##  [863,] 1.0084409 0.9502844
    ##  [864,] 0.9616466 1.0056533
    ##  [865,] 0.9959834 1.0810569
    ##  [866,] 0.9683666 0.9824692
    ##  [867,] 0.8987527 1.0641582
    ##  [868,] 0.9160065 1.0473785
    ##  [869,] 0.8717188 0.9697657
    ##  [870,] 1.0931577 0.9667451
    ##  [871,] 1.2051432 1.0197571
    ##  [872,] 1.0866577 0.9664647
    ##  [873,] 0.8672311 1.1174151
    ##  [874,] 0.8800225 0.9877608
    ##  [875,] 0.9443667 1.2427277
    ##  [876,] 0.8939368 0.9375530
    ##  [877,] 1.0055776 0.9814461
    ##  [878,] 0.9447647 1.0967718
    ##  [879,] 1.0230638 0.9452838
    ##  [880,] 0.9299225 1.0707332
    ##  [881,] 1.0049603 1.1095836
    ##  [882,] 0.9764401 0.9234940
    ##  [883,] 0.9811320 1.0477180
    ##  [884,] 0.8860742 0.9832692
    ##  [885,] 0.9921751 0.9608090
    ##  [886,] 1.0435646 1.0600318
    ##  [887,] 1.1106368 1.0429525
    ##  [888,] 0.9029847 1.0747693
    ##  [889,] 1.0517654 1.1043093
    ##  [890,] 0.9300660 1.0335110
    ##  [891,] 0.9447687 1.0108163
    ##  [892,] 0.8641682 1.0666088
    ##  [893,] 0.9213797 1.0887232
    ##  [894,] 0.9827931 1.0527900
    ##  [895,] 1.1463204 1.0464140
    ##  [896,] 1.0242059 0.8968239
    ##  [897,] 1.0358074 0.9323151
    ##  [898,] 0.8098828 1.0201859
    ##  [899,] 0.9660818 1.0453540
    ##  [900,] 0.9257626 0.9898941
    ##  [901,] 0.9406827 1.0394911
    ##  [902,] 0.9802793 1.0387884
    ##  [903,] 1.1982388 0.9761859
    ##  [904,] 1.0527154 0.8936029
    ##  [905,] 0.9253210 1.0162115
    ##  [906,] 0.9537891 0.9855714
    ##  [907,] 1.1358189 0.8917485
    ##  [908,] 0.9615261 1.1244603
    ##  [909,] 1.0329919 0.9160823
    ##  [910,] 1.0382319 0.9151813
    ##  [911,] 1.0521275 0.9849822
    ##  [912,] 0.9152559 1.0505183
    ##  [913,] 0.9349982 0.9456630
    ##  [914,] 0.9418645 1.0397835
    ##  [915,] 1.0257449 1.0127150
    ##  [916,] 0.9192848 1.0324150
    ##  [917,] 0.8903530 0.9580320
    ##  [918,] 1.0372470 1.0532542
    ##  [919,] 1.0789425 1.1166084
    ##  [920,] 0.8633820 0.8517114
    ##  [921,] 0.8939772 0.9660395
    ##  [922,] 1.0289966 1.0914506
    ##  [923,] 1.0191781 1.0670840
    ##  [924,] 0.9534113 1.0188594
    ##  [925,] 0.9639898 1.0889665
    ##  [926,] 0.9932369 0.9865375
    ##  [927,] 0.9914002 0.8993844
    ##  [928,] 1.0953473 0.9776709
    ##  [929,] 1.0243304 0.9634915
    ##  [930,] 1.0737289 0.9980714
    ##  [931,] 0.8493196 0.8776111
    ##  [932,] 1.0337207 0.8808628
    ##  [933,] 0.9077633 1.0132058
    ##  [934,] 1.0401527 1.0788335
    ##  [935,] 0.9728057 0.9190930
    ##  [936,] 0.9136534 0.9681767
    ##  [937,] 0.9718410 0.9734074
    ##  [938,] 0.8013217 0.9551793
    ##  [939,] 1.2576565 1.0573464
    ##  [940,] 0.9311768 1.0725939
    ##  [941,] 0.9487521 0.9874652
    ##  [942,] 1.3560823 0.9977878
    ##  [943,] 0.8407918 0.9224680
    ##  [944,] 1.1373721 0.8994977
    ##  [945,] 0.9143340 0.9515472
    ##  [946,] 0.9836673 1.0582391
    ##  [947,] 0.9926987 0.9717312
    ##  [948,] 1.2058899 1.0170934
    ##  [949,] 1.0578684 0.9046883
    ##  [950,] 0.8424216 0.9399182
    ##  [951,] 0.9171803 1.0895950
    ##  [952,] 1.0804320 0.9875105
    ##  [953,] 1.0457073 0.9375547
    ##  [954,] 0.9965746 0.9323030
    ##  [955,] 1.0244696 1.0180279
    ##  [956,] 1.0409610 1.0252629
    ##  [957,] 0.9113145 0.9652427
    ##  [958,] 1.1666450 1.0184796
    ##  [959,] 0.8755509 0.8775956
    ##  [960,] 0.8869760 0.9318794
    ##  [961,] 1.0758201 1.0099227
    ##  [962,] 0.9208705 0.9934749
    ##  [963,] 1.0249357 0.9311806
    ##  [964,] 1.0710522 1.0132168
    ##  [965,] 0.8749687 0.9746530
    ##  [966,] 0.9859856 0.9699251
    ##  [967,] 1.1347164 0.9601095
    ##  [968,] 0.9572813 1.0440613
    ##  [969,] 0.9150124 0.9674912
    ##  [970,] 0.9451598 0.9731706
    ##  [971,] 0.8742447 1.0711193
    ##  [972,] 0.8326364 1.1185858
    ##  [973,] 1.0351203 0.9643120
    ##  [974,] 1.0753120 0.9210882
    ##  [975,] 0.9079309 1.0909104
    ##  [976,] 1.0330600 1.0150539
    ##  [977,] 0.8856314 0.8955992
    ##  [978,] 1.2640872 1.0889025
    ##  [979,] 1.1483564 0.9931464
    ##  [980,] 0.8970450 0.9779669
    ##  [981,] 0.9017425 0.9457510
    ##  [982,] 0.8939568 0.9403775
    ##  [983,] 0.8867260 1.1286917
    ##  [984,] 1.0058375 1.1432302
    ##  [985,] 1.0964319 0.9728180
    ##  [986,] 1.1091314 0.9283088
    ##  [987,] 0.9991410 0.9296889
    ##  [988,] 1.0627876 0.8973085
    ##  [989,] 0.9869603 0.9454768
    ##  [990,] 1.0098333 1.1313193
    ##  [991,] 0.8705572 1.0416625
    ##  [992,] 0.9253594 0.9500899
    ##  [993,] 0.8883746 0.9853059
    ##  [994,] 0.9895588 1.0354198
    ##  [995,] 0.8887422 1.0271998
    ##  [996,] 0.8816639 1.0258848
    ##  [997,] 1.0013238 1.1710159
    ##  [998,] 0.9176645 1.0492721
    ##  [999,] 0.9046724 0.9725742
    ## 
    ## $model.matrix
    ##    (Intercept) site1 pH1 pH2
    ## 1            1     1   1   0
    ## 2            1     1   1   0
    ## 3            1     1   1   0
    ## 4            1     1   1   0
    ## 5            1     1  -1  -1
    ## 6            1     1  -1  -1
    ## 7            1     1  -1  -1
    ## 8            1     1  -1  -1
    ## 9            1     1   0   1
    ## 10           1     1   0   1
    ## 11           1     1   0   1
    ## 12           1     1   0   1
    ## 13           1    -1   1   0
    ## 14           1    -1   1   0
    ## 15           1    -1   1   0
    ## 16           1    -1   1   0
    ## 17           1    -1  -1  -1
    ## 18           1    -1  -1  -1
    ## 19           1    -1  -1  -1
    ## 20           1    -1  -1  -1
    ## 
    ## $terms
    ## dist_tab_assay ~ site + pH
    ## attr(,"variables")
    ## list(dist_tab_assay, site, pH)
    ## attr(,"factors")
    ##                site pH
    ## dist_tab_assay    0  0
    ## site              1  0
    ## pH                0  1
    ## attr(,"term.labels")
    ## [1] "site" "pH"  
    ## attr(,"order")
    ## [1] 1 1
    ## attr(,"intercept")
    ## [1] 1
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>
    ## 
    ## attr(,"class")
    ## [1] "adonis"

``` r
anova(betadisper(dist_tab_assay,samples$site))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df  Sum Sq Mean Sq F value  Pr(>F)  
    ## Groups     1  402.36  402.36  3.2462 0.08836 .
    ## Residuals 18 2231.02  123.95                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(betadisper(dist_tab_assay,samples$pH))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value   Pr(>F)   
    ## Groups     2 1219.4  609.70  6.6626 0.007303 **
    ## Residuals 17 1555.7   91.51                    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
count_tab_assay <- assay(vsdNatSim)
dist_tab_assay <- dist(t(count_tab_assay),method="euclidian")
adonis(data=samplesNatSim,dist_tab_assay ~ site_pH, method="euclidian")
```

    ## 'adonis' will be deprecated: use 'adonis2' instead

    ## $aov.tab
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## Terms added sequentially (first to last)
    ## 
    ##           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
    ## site_pH    1     28338   28338  1.5518 0.20549  0.026 *
    ## Residuals  6    109566   18261         0.79451         
    ## Total      7    137904                 1.00000         
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## $call
    ## adonis(formula = dist_tab_assay ~ site_pH, data = samplesNatSim, 
    ##     method = "euclidian")
    ## 
    ## $coefficients
    ## NULL
    ## 
    ## $coef.sites
    ##                  [,1]     [,2]      [,3]      [,4]     [,5]      [,6]      [,7]
    ## (Intercept) 174.40097 188.2578 174.91784 170.54814 173.6810 167.39969 166.81613
    ## site_pH1    -29.13536 -33.0374 -22.50114 -21.23729  35.3239  29.84631  35.52947
    ##                  [,8]
    ## (Intercept) 169.89554
    ## site_pH1     35.54383
    ## 
    ## $f.perms
    ##              [,1]
    ##    [1,] 0.9715200
    ##    [2,] 0.8499509
    ##    [3,] 1.0014790
    ##    [4,] 1.1486690
    ##    [5,] 0.9390040
    ##    [6,] 1.2768584
    ##    [7,] 1.0735319
    ##    [8,] 0.9121887
    ##    [9,] 1.5518528
    ##   [10,] 0.8949477
    ##   [11,] 0.8386694
    ##   [12,] 0.8313592
    ##   [13,] 1.2768584
    ##   [14,] 0.9121887
    ##   [15,] 0.8386694
    ##   [16,] 1.0008314
    ##   [17,] 0.9349355
    ##   [18,] 0.8759996
    ##   [19,] 0.9107817
    ##   [20,] 0.8595390
    ##   [21,] 1.1486690
    ##   [22,] 1.1372653
    ##   [23,] 1.1233072
    ##   [24,] 1.3379525
    ##   [25,] 0.8759996
    ##   [26,] 0.8646186
    ##   [27,] 1.0748758
    ##   [28,] 0.8856563
    ##   [29,] 0.8001891
    ##   [30,] 0.8313592
    ##   [31,] 1.2768584
    ##   [32,] 1.0498496
    ##   [33,] 0.9349355
    ##   [34,] 0.8080195
    ##   [35,] 0.8499509
    ##   [36,] 0.9107817
    ##   [37,] 0.9390040
    ##   [38,] 1.0498496
    ##   [39,] 1.0498496
    ##   [40,] 1.1372653
    ##   [41,] 1.0498496
    ##   [42,] 0.9593707
    ##   [43,] 1.0008314
    ##   [44,] 1.0014790
    ##   [45,] 0.8313592
    ##   [46,] 0.9390040
    ##   [47,] 1.1002390
    ##   [48,] 0.9593707
    ##   [49,] 0.9107817
    ##   [50,] 1.0252180
    ##   [51,] 1.0735319
    ##   [52,] 1.0498496
    ##   [53,] 1.3379525
    ##   [54,] 1.0498496
    ##   [55,] 0.8080195
    ##   [56,] 1.5518528
    ##   [57,] 1.1002390
    ##   [58,] 0.8646186
    ##   [59,] 0.8499509
    ##   [60,] 0.8759996
    ##   [61,] 0.9715200
    ##   [62,] 1.0498496
    ##   [63,] 1.0440420
    ##   [64,] 0.8386694
    ##   [65,] 0.9715200
    ##   [66,] 0.8595390
    ##   [67,] 0.9715200
    ##   [68,] 0.8080195
    ##   [69,] 0.8759996
    ##   [70,] 1.0252180
    ##   [71,] 0.8499509
    ##   [72,] 1.0440420
    ##   [73,] 0.8595390
    ##   [74,] 0.8804854
    ##   [75,] 1.3379525
    ##   [76,] 1.0008314
    ##   [77,] 1.3379525
    ##   [78,] 1.3379525
    ##   [79,] 0.8856563
    ##   [80,] 0.8001891
    ##   [81,] 0.8968552
    ##   [82,] 0.9715200
    ##   [83,] 0.8968552
    ##   [84,] 0.8804854
    ##   [85,] 1.1233072
    ##   [86,] 0.8313592
    ##   [87,] 1.0252180
    ##   [88,] 0.8080195
    ##   [89,] 0.9715200
    ##   [90,] 0.9593707
    ##   [91,] 0.9107817
    ##   [92,] 1.0735319
    ##   [93,] 0.9121887
    ##   [94,] 1.0252180
    ##   [95,] 0.8313592
    ##   [96,] 0.8080195
    ##   [97,] 1.0498496
    ##   [98,] 1.0252180
    ##   [99,] 0.8080195
    ##  [100,] 1.1641895
    ##  [101,] 0.8499509
    ##  [102,] 1.3379525
    ##  [103,] 0.8646186
    ##  [104,] 0.9107817
    ##  [105,] 1.1233072
    ##  [106,] 0.8949477
    ##  [107,] 1.0440420
    ##  [108,] 0.9121887
    ##  [109,] 1.3379525
    ##  [110,] 1.1372653
    ##  [111,] 0.8001891
    ##  [112,] 0.8949477
    ##  [113,] 1.1641895
    ##  [114,] 1.0735319
    ##  [115,] 0.8499509
    ##  [116,] 0.8856563
    ##  [117,] 0.9107817
    ##  [118,] 0.8080195
    ##  [119,] 1.2768584
    ##  [120,] 1.0748758
    ##  [121,] 0.8313592
    ##  [122,] 1.0498496
    ##  [123,] 1.0008314
    ##  [124,] 1.1372653
    ##  [125,] 1.5518528
    ##  [126,] 1.2167857
    ##  [127,] 0.8949477
    ##  [128,] 1.1002390
    ##  [129,] 0.8949477
    ##  [130,] 0.8595390
    ##  [131,] 0.8759996
    ##  [132,] 1.0014790
    ##  [133,] 1.0014790
    ##  [134,] 0.8646186
    ##  [135,] 1.2167857
    ##  [136,] 1.2768584
    ##  [137,] 0.9107817
    ##  [138,] 1.2768584
    ##  [139,] 1.1002390
    ##  [140,] 0.8856563
    ##  [141,] 1.5518528
    ##  [142,] 1.2768584
    ##  [143,] 1.2768584
    ##  [144,] 1.0014790
    ##  [145,] 0.8949477
    ##  [146,] 0.9715200
    ##  [147,] 1.1233072
    ##  [148,] 0.8313592
    ##  [149,] 0.9349355
    ##  [150,] 0.9390040
    ##  [151,] 0.8080195
    ##  [152,] 1.2167857
    ##  [153,] 1.0748758
    ##  [154,] 0.8001891
    ##  [155,] 0.8001891
    ##  [156,] 1.1641895
    ##  [157,] 1.1372653
    ##  [158,] 0.9715200
    ##  [159,] 0.8759996
    ##  [160,] 1.1486690
    ##  [161,] 1.1372653
    ##  [162,] 0.9121887
    ##  [163,] 1.1486690
    ##  [164,] 0.8001891
    ##  [165,] 0.8499509
    ##  [166,] 0.8949477
    ##  [167,] 0.8499509
    ##  [168,] 0.8901779
    ##  [169,] 1.0008314
    ##  [170,] 0.8386694
    ##  [171,] 0.8001891
    ##  [172,] 0.8080195
    ##  [173,] 0.8968552
    ##  [174,] 1.0735319
    ##  [175,] 0.8499509
    ##  [176,] 0.8080195
    ##  [177,] 0.8499509
    ##  [178,] 1.0748758
    ##  [179,] 1.0014790
    ##  [180,] 1.0440420
    ##  [181,] 0.8386694
    ##  [182,] 0.9349355
    ##  [183,] 1.0735319
    ##  [184,] 0.8499509
    ##  [185,] 0.9390040
    ##  [186,] 0.8499509
    ##  [187,] 0.8001891
    ##  [188,] 1.3379525
    ##  [189,] 0.8856563
    ##  [190,] 0.8759996
    ##  [191,] 0.8595390
    ##  [192,] 1.3379525
    ##  [193,] 1.5518528
    ##  [194,] 0.8759996
    ##  [195,] 0.8313592
    ##  [196,] 0.8949477
    ##  [197,] 0.9390040
    ##  [198,] 0.9121887
    ##  [199,] 0.8646186
    ##  [200,] 0.9715200
    ##  [201,] 0.8080195
    ##  [202,] 0.9121887
    ##  [203,] 0.8804854
    ##  [204,] 1.0252180
    ##  [205,] 1.0014790
    ##  [206,] 1.0735319
    ##  [207,] 1.0748758
    ##  [208,] 1.5518528
    ##  [209,] 1.3379525
    ##  [210,] 0.8949477
    ##  [211,] 0.8804854
    ##  [212,] 0.8804854
    ##  [213,] 1.0014790
    ##  [214,] 0.9349355
    ##  [215,] 0.9121887
    ##  [216,] 0.9121887
    ##  [217,] 0.9349355
    ##  [218,] 0.8595390
    ##  [219,] 0.8499509
    ##  [220,] 0.8001891
    ##  [221,] 0.9349355
    ##  [222,] 1.0498496
    ##  [223,] 1.0252180
    ##  [224,] 0.8001891
    ##  [225,] 1.1641895
    ##  [226,] 1.0014790
    ##  [227,] 0.9349355
    ##  [228,] 1.5518528
    ##  [229,] 0.8901779
    ##  [230,] 0.8759996
    ##  [231,] 1.1233072
    ##  [232,] 0.8080195
    ##  [233,] 0.8968552
    ##  [234,] 1.1002390
    ##  [235,] 0.8313592
    ##  [236,] 1.1372653
    ##  [237,] 1.0748758
    ##  [238,] 0.8804854
    ##  [239,] 0.8080195
    ##  [240,] 0.9715200
    ##  [241,] 1.1641895
    ##  [242,] 1.2768584
    ##  [243,] 0.8313592
    ##  [244,] 0.8804854
    ##  [245,] 1.1486690
    ##  [246,] 0.8949477
    ##  [247,] 0.9390040
    ##  [248,] 0.8949477
    ##  [249,] 1.0014790
    ##  [250,] 0.8386694
    ##  [251,] 0.8595390
    ##  [252,] 1.0735319
    ##  [253,] 0.9593707
    ##  [254,] 1.5518528
    ##  [255,] 0.8001891
    ##  [256,] 0.8804854
    ##  [257,] 1.1372653
    ##  [258,] 1.0252180
    ##  [259,] 1.0008314
    ##  [260,] 0.8901779
    ##  [261,] 1.0498496
    ##  [262,] 1.2768584
    ##  [263,] 1.0735319
    ##  [264,] 1.2167857
    ##  [265,] 1.0014790
    ##  [266,] 0.9390040
    ##  [267,] 1.1372653
    ##  [268,] 1.0735319
    ##  [269,] 0.8595390
    ##  [270,] 1.3379525
    ##  [271,] 0.8499509
    ##  [272,] 1.0735319
    ##  [273,] 1.3379525
    ##  [274,] 1.2167857
    ##  [275,] 0.8759996
    ##  [276,] 0.8313592
    ##  [277,] 0.8901779
    ##  [278,] 0.8386694
    ##  [279,] 1.1486690
    ##  [280,] 0.9593707
    ##  [281,] 0.8856563
    ##  [282,] 1.1233072
    ##  [283,] 0.9349355
    ##  [284,] 1.1233072
    ##  [285,] 1.1233072
    ##  [286,] 0.8313592
    ##  [287,] 0.8646186
    ##  [288,] 1.5518528
    ##  [289,] 1.0014790
    ##  [290,] 0.8499509
    ##  [291,] 1.1233072
    ##  [292,] 1.2167857
    ##  [293,] 0.9390040
    ##  [294,] 0.8386694
    ##  [295,] 0.8804854
    ##  [296,] 1.0735319
    ##  [297,] 1.1641895
    ##  [298,] 0.8646186
    ##  [299,] 1.1641895
    ##  [300,] 1.3379525
    ##  [301,] 0.8499509
    ##  [302,] 0.9121887
    ##  [303,] 1.0735319
    ##  [304,] 0.9121887
    ##  [305,] 0.8499509
    ##  [306,] 0.8804854
    ##  [307,] 1.1372653
    ##  [308,] 0.9121887
    ##  [309,] 0.8313592
    ##  [310,] 1.3379525
    ##  [311,] 0.8001891
    ##  [312,] 0.9121887
    ##  [313,] 0.8901779
    ##  [314,] 1.2167857
    ##  [315,] 0.9390040
    ##  [316,] 0.8856563
    ##  [317,] 1.1641895
    ##  [318,] 1.0748758
    ##  [319,] 0.8949477
    ##  [320,] 1.0748758
    ##  [321,] 1.0498496
    ##  [322,] 0.9593707
    ##  [323,] 0.8646186
    ##  [324,] 0.9349355
    ##  [325,] 0.8804854
    ##  [326,] 0.8001891
    ##  [327,] 0.9107817
    ##  [328,] 0.8386694
    ##  [329,] 1.1486690
    ##  [330,] 0.8968552
    ##  [331,] 1.0498496
    ##  [332,] 1.1233072
    ##  [333,] 1.1486690
    ##  [334,] 1.0008314
    ##  [335,] 0.8901779
    ##  [336,] 0.8759996
    ##  [337,] 1.1486690
    ##  [338,] 0.8901779
    ##  [339,] 0.8001891
    ##  [340,] 0.9715200
    ##  [341,] 1.0440420
    ##  [342,] 1.0440420
    ##  [343,] 1.5518528
    ##  [344,] 0.8001891
    ##  [345,] 1.1486690
    ##  [346,] 1.1641895
    ##  [347,] 0.9593707
    ##  [348,] 0.8001891
    ##  [349,] 1.1486690
    ##  [350,] 1.2768584
    ##  [351,] 1.1486690
    ##  [352,] 1.0252180
    ##  [353,] 0.8001891
    ##  [354,] 0.9121887
    ##  [355,] 0.8499509
    ##  [356,] 0.8968552
    ##  [357,] 0.8313592
    ##  [358,] 0.8499509
    ##  [359,] 1.1372653
    ##  [360,] 1.0014790
    ##  [361,] 1.2167857
    ##  [362,] 0.8901779
    ##  [363,] 0.9349355
    ##  [364,] 1.2768584
    ##  [365,] 1.1233072
    ##  [366,] 1.0252180
    ##  [367,] 0.8949477
    ##  [368,] 0.8804854
    ##  [369,] 1.1002390
    ##  [370,] 0.8949477
    ##  [371,] 0.9107817
    ##  [372,] 0.8646186
    ##  [373,] 0.9593707
    ##  [374,] 1.1372653
    ##  [375,] 0.8595390
    ##  [376,] 1.1233072
    ##  [377,] 0.8313592
    ##  [378,] 0.8759996
    ##  [379,] 0.8901779
    ##  [380,] 1.3379525
    ##  [381,] 1.0748758
    ##  [382,] 1.0252180
    ##  [383,] 0.9107817
    ##  [384,] 0.8499509
    ##  [385,] 1.2167857
    ##  [386,] 0.8386694
    ##  [387,] 0.8968552
    ##  [388,] 0.9593707
    ##  [389,] 1.2167857
    ##  [390,] 0.9715200
    ##  [391,] 0.9715200
    ##  [392,] 1.2768584
    ##  [393,] 1.0748758
    ##  [394,] 0.9121887
    ##  [395,] 0.8386694
    ##  [396,] 0.8901779
    ##  [397,] 1.1372653
    ##  [398,] 1.1641895
    ##  [399,] 0.8968552
    ##  [400,] 1.1372653
    ##  [401,] 1.0252180
    ##  [402,] 1.3379525
    ##  [403,] 0.8499509
    ##  [404,] 0.8949477
    ##  [405,] 1.1233072
    ##  [406,] 1.0440420
    ##  [407,] 1.1002390
    ##  [408,] 0.8595390
    ##  [409,] 1.2768584
    ##  [410,] 1.1002390
    ##  [411,] 1.0498496
    ##  [412,] 1.1641895
    ##  [413,] 0.8499509
    ##  [414,] 1.5518528
    ##  [415,] 0.8386694
    ##  [416,] 1.3379525
    ##  [417,] 1.0014790
    ##  [418,] 0.8901779
    ##  [419,] 0.9121887
    ##  [420,] 1.0440420
    ##  [421,] 0.9593707
    ##  [422,] 0.8499509
    ##  [423,] 1.0252180
    ##  [424,] 0.8080195
    ##  [425,] 0.8595390
    ##  [426,] 1.1486690
    ##  [427,] 1.0014790
    ##  [428,] 0.8759996
    ##  [429,] 0.8968552
    ##  [430,] 0.8313592
    ##  [431,] 0.8646186
    ##  [432,] 1.3379525
    ##  [433,] 0.9390040
    ##  [434,] 0.8499509
    ##  [435,] 0.8313592
    ##  [436,] 0.8646186
    ##  [437,] 1.2768584
    ##  [438,] 0.8804854
    ##  [439,] 0.9593707
    ##  [440,] 1.0008314
    ##  [441,] 0.8080195
    ##  [442,] 1.3379525
    ##  [443,] 1.5518528
    ##  [444,] 0.8646186
    ##  [445,] 1.2167857
    ##  [446,] 0.8646186
    ##  [447,] 1.1372653
    ##  [448,] 0.8759996
    ##  [449,] 1.1233072
    ##  [450,] 1.0748758
    ##  [451,] 1.0014790
    ##  [452,] 1.0748758
    ##  [453,] 1.1486690
    ##  [454,] 1.1641895
    ##  [455,] 0.8804854
    ##  [456,] 0.8949477
    ##  [457,] 0.8080195
    ##  [458,] 0.8313592
    ##  [459,] 1.0748758
    ##  [460,] 0.8499509
    ##  [461,] 1.0252180
    ##  [462,] 1.0252180
    ##  [463,] 0.9593707
    ##  [464,] 0.8313592
    ##  [465,] 0.8499509
    ##  [466,] 0.8001891
    ##  [467,] 1.1641895
    ##  [468,] 0.8759996
    ##  [469,] 1.0735319
    ##  [470,] 0.8856563
    ##  [471,] 1.0440420
    ##  [472,] 1.1372653
    ##  [473,] 1.2768584
    ##  [474,] 1.2768584
    ##  [475,] 0.8968552
    ##  [476,] 0.8001891
    ##  [477,] 1.0748758
    ##  [478,] 1.2768584
    ##  [479,] 0.9349355
    ##  [480,] 0.9390040
    ##  [481,] 0.8968552
    ##  [482,] 0.8386694
    ##  [483,] 0.8968552
    ##  [484,] 0.8499509
    ##  [485,] 1.1233072
    ##  [486,] 1.0440420
    ##  [487,] 0.8595390
    ##  [488,] 1.1002390
    ##  [489,] 0.8804854
    ##  [490,] 0.8949477
    ##  [491,] 0.9715200
    ##  [492,] 1.5518528
    ##  [493,] 1.1233072
    ##  [494,] 0.8804854
    ##  [495,] 1.2167857
    ##  [496,] 0.8901779
    ##  [497,] 1.0252180
    ##  [498,] 1.1233072
    ##  [499,] 1.0014790
    ##  [500,] 0.8804854
    ##  [501,] 1.1641895
    ##  [502,] 1.3379525
    ##  [503,] 1.1372653
    ##  [504,] 0.8313592
    ##  [505,] 0.8856563
    ##  [506,] 0.8001891
    ##  [507,] 0.8001891
    ##  [508,] 0.8901779
    ##  [509,] 1.1233072
    ##  [510,] 0.9121887
    ##  [511,] 0.8759996
    ##  [512,] 1.2768584
    ##  [513,] 0.8949477
    ##  [514,] 1.0440420
    ##  [515,] 0.8759996
    ##  [516,] 0.8968552
    ##  [517,] 1.1641895
    ##  [518,] 1.0014790
    ##  [519,] 0.8949477
    ##  [520,] 0.8080195
    ##  [521,] 0.8080195
    ##  [522,] 0.8080195
    ##  [523,] 0.8386694
    ##  [524,] 0.8949477
    ##  [525,] 0.8646186
    ##  [526,] 0.8646186
    ##  [527,] 0.8595390
    ##  [528,] 1.0008314
    ##  [529,] 0.8759996
    ##  [530,] 1.0440420
    ##  [531,] 1.0252180
    ##  [532,] 1.0498496
    ##  [533,] 1.0440420
    ##  [534,] 0.9715200
    ##  [535,] 0.8499509
    ##  [536,] 0.9107817
    ##  [537,] 0.9390040
    ##  [538,] 0.9593707
    ##  [539,] 0.8499509
    ##  [540,] 1.1641895
    ##  [541,] 1.2768584
    ##  [542,] 0.8646186
    ##  [543,] 1.0440420
    ##  [544,] 0.8646186
    ##  [545,] 1.0440420
    ##  [546,] 0.8499509
    ##  [547,] 0.8968552
    ##  [548,] 0.8080195
    ##  [549,] 1.3379525
    ##  [550,] 1.3379525
    ##  [551,] 0.8001891
    ##  [552,] 1.1002390
    ##  [553,] 0.8968552
    ##  [554,] 0.9593707
    ##  [555,] 0.8759996
    ##  [556,] 0.8499509
    ##  [557,] 1.0014790
    ##  [558,] 0.8804854
    ##  [559,] 1.1641895
    ##  [560,] 1.0252180
    ##  [561,] 0.9121887
    ##  [562,] 0.8499509
    ##  [563,] 1.0252180
    ##  [564,] 0.8646186
    ##  [565,] 1.0735319
    ##  [566,] 0.9390040
    ##  [567,] 1.2167857
    ##  [568,] 0.8759996
    ##  [569,] 0.9593707
    ##  [570,] 0.9593707
    ##  [571,] 1.3379525
    ##  [572,] 1.1641895
    ##  [573,] 0.8968552
    ##  [574,] 1.1233072
    ##  [575,] 0.9121887
    ##  [576,] 1.0748758
    ##  [577,] 1.1002390
    ##  [578,] 0.9715200
    ##  [579,] 1.1486690
    ##  [580,] 1.1233072
    ##  [581,] 0.8968552
    ##  [582,] 0.9349355
    ##  [583,] 0.8759996
    ##  [584,] 1.3379525
    ##  [585,] 0.8804854
    ##  [586,] 1.1641895
    ##  [587,] 0.9121887
    ##  [588,] 1.1002390
    ##  [589,] 0.8595390
    ##  [590,] 0.9715200
    ##  [591,] 0.8001891
    ##  [592,] 0.8001891
    ##  [593,] 0.8313592
    ##  [594,] 1.0748758
    ##  [595,] 0.9107817
    ##  [596,] 1.5518528
    ##  [597,] 0.8804854
    ##  [598,] 0.9593707
    ##  [599,] 1.0498496
    ##  [600,] 1.0735319
    ##  [601,] 0.9349355
    ##  [602,] 0.9121887
    ##  [603,] 0.8804854
    ##  [604,] 1.1233072
    ##  [605,] 1.0748758
    ##  [606,] 0.8804854
    ##  [607,] 0.8080195
    ##  [608,] 0.8499509
    ##  [609,] 0.8001891
    ##  [610,] 0.9121887
    ##  [611,] 1.5518528
    ##  [612,] 0.8646186
    ##  [613,] 0.8759996
    ##  [614,] 0.8759996
    ##  [615,] 1.2167857
    ##  [616,] 1.0008314
    ##  [617,] 0.9390040
    ##  [618,] 0.8001891
    ##  [619,] 1.0008314
    ##  [620,] 0.8499509
    ##  [621,] 0.8313592
    ##  [622,] 0.8856563
    ##  [623,] 0.8595390
    ##  [624,] 1.0008314
    ##  [625,] 1.0735319
    ##  [626,] 1.1233072
    ##  [627,] 0.9349355
    ##  [628,] 1.0735319
    ##  [629,] 0.8949477
    ##  [630,] 1.3379525
    ##  [631,] 0.8804854
    ##  [632,] 0.8901779
    ##  [633,] 0.9349355
    ##  [634,] 0.9349355
    ##  [635,] 0.9121887
    ##  [636,] 0.8949477
    ##  [637,] 1.0735319
    ##  [638,] 0.9121887
    ##  [639,] 1.1486690
    ##  [640,] 0.8595390
    ##  [641,] 1.0252180
    ##  [642,] 1.0440420
    ##  [643,] 1.3379525
    ##  [644,] 1.1372653
    ##  [645,] 1.5518528
    ##  [646,] 1.0014790
    ##  [647,] 1.3379525
    ##  [648,] 0.8804854
    ##  [649,] 0.8804854
    ##  [650,] 0.8901779
    ##  [651,] 0.9715200
    ##  [652,] 1.0014790
    ##  [653,] 0.8901779
    ##  [654,] 1.1233072
    ##  [655,] 1.2768584
    ##  [656,] 1.0008314
    ##  [657,] 0.8313592
    ##  [658,] 0.8646186
    ##  [659,] 1.3379525
    ##  [660,] 1.0014790
    ##  [661,] 1.0735319
    ##  [662,] 0.9121887
    ##  [663,] 0.8646186
    ##  [664,] 1.0735319
    ##  [665,] 0.9107817
    ##  [666,] 1.0014790
    ##  [667,] 1.2768584
    ##  [668,] 0.9715200
    ##  [669,] 0.8595390
    ##  [670,] 0.8499509
    ##  [671,] 0.9715200
    ##  [672,] 1.0748758
    ##  [673,] 0.8759996
    ##  [674,] 0.8499509
    ##  [675,] 0.9593707
    ##  [676,] 0.8313592
    ##  [677,] 1.0498496
    ##  [678,] 1.2167857
    ##  [679,] 0.9121887
    ##  [680,] 1.1372653
    ##  [681,] 1.2167857
    ##  [682,] 0.9390040
    ##  [683,] 0.8646186
    ##  [684,] 1.2768584
    ##  [685,] 0.8595390
    ##  [686,] 1.5518528
    ##  [687,] 0.9107817
    ##  [688,] 1.3379525
    ##  [689,] 1.5518528
    ##  [690,] 0.9107817
    ##  [691,] 1.1372653
    ##  [692,] 1.0748758
    ##  [693,] 0.8001891
    ##  [694,] 1.1002390
    ##  [695,] 0.8949477
    ##  [696,] 0.8949477
    ##  [697,] 1.2768584
    ##  [698,] 0.9121887
    ##  [699,] 0.8968552
    ##  [700,] 0.8949477
    ##  [701,] 1.5518528
    ##  [702,] 1.0014790
    ##  [703,] 1.0008314
    ##  [704,] 0.8856563
    ##  [705,] 0.8901779
    ##  [706,] 0.9121887
    ##  [707,] 0.8499509
    ##  [708,] 0.9715200
    ##  [709,] 0.8386694
    ##  [710,] 0.8968552
    ##  [711,] 1.2167857
    ##  [712,] 1.2167857
    ##  [713,] 0.8595390
    ##  [714,] 1.1486690
    ##  [715,] 1.0748758
    ##  [716,] 0.9390040
    ##  [717,] 0.9593707
    ##  [718,] 1.0735319
    ##  [719,] 0.9390040
    ##  [720,] 1.1486690
    ##  [721,] 0.8804854
    ##  [722,] 0.8595390
    ##  [723,] 0.9593707
    ##  [724,] 1.2167857
    ##  [725,] 0.9349355
    ##  [726,] 0.8386694
    ##  [727,] 1.2768584
    ##  [728,] 1.0252180
    ##  [729,] 0.8759996
    ##  [730,] 0.8386694
    ##  [731,] 1.0008314
    ##  [732,] 0.9715200
    ##  [733,] 0.9349355
    ##  [734,] 0.8001891
    ##  [735,] 1.3379525
    ##  [736,] 0.9121887
    ##  [737,] 1.2167857
    ##  [738,] 1.1641895
    ##  [739,] 0.9107817
    ##  [740,] 0.8646186
    ##  [741,] 1.0498496
    ##  [742,] 0.8313592
    ##  [743,] 0.8595390
    ##  [744,] 1.0748758
    ##  [745,] 0.8804854
    ##  [746,] 0.8595390
    ##  [747,] 0.8595390
    ##  [748,] 0.9107817
    ##  [749,] 0.8001891
    ##  [750,] 0.8001891
    ##  [751,] 0.9349355
    ##  [752,] 1.1641895
    ##  [753,] 0.8386694
    ##  [754,] 0.8804854
    ##  [755,] 0.8313592
    ##  [756,] 1.2768584
    ##  [757,] 1.0735319
    ##  [758,] 1.1641895
    ##  [759,] 1.1486690
    ##  [760,] 0.9121887
    ##  [761,] 0.8313592
    ##  [762,] 0.9349355
    ##  [763,] 1.1372653
    ##  [764,] 0.8901779
    ##  [765,] 0.8001891
    ##  [766,] 0.8759996
    ##  [767,] 0.8001891
    ##  [768,] 1.0014790
    ##  [769,] 1.0008314
    ##  [770,] 0.9715200
    ##  [771,] 0.8313592
    ##  [772,] 1.5518528
    ##  [773,] 1.1641895
    ##  [774,] 1.1233072
    ##  [775,] 0.8856563
    ##  [776,] 0.8386694
    ##  [777,] 1.1641895
    ##  [778,] 0.8646186
    ##  [779,] 1.1233072
    ##  [780,] 1.0008314
    ##  [781,] 1.2167857
    ##  [782,] 1.0014790
    ##  [783,] 1.0498496
    ##  [784,] 1.1641895
    ##  [785,] 1.1486690
    ##  [786,] 1.0014790
    ##  [787,] 0.8386694
    ##  [788,] 0.8968552
    ##  [789,] 0.8968552
    ##  [790,] 1.0498496
    ##  [791,] 0.8856563
    ##  [792,] 0.8759996
    ##  [793,] 1.3379525
    ##  [794,] 1.2768584
    ##  [795,] 1.0748758
    ##  [796,] 1.0498496
    ##  [797,] 1.0748758
    ##  [798,] 0.8968552
    ##  [799,] 0.8804854
    ##  [800,] 0.8901779
    ##  [801,] 1.0748758
    ##  [802,] 1.2768584
    ##  [803,] 1.0014790
    ##  [804,] 1.5518528
    ##  [805,] 1.0014790
    ##  [806,] 0.9349355
    ##  [807,] 0.8595390
    ##  [808,] 0.8856563
    ##  [809,] 0.8313592
    ##  [810,] 0.9715200
    ##  [811,] 0.8386694
    ##  [812,] 0.9121887
    ##  [813,] 0.8080195
    ##  [814,] 1.1372653
    ##  [815,] 0.8386694
    ##  [816,] 0.8968552
    ##  [817,] 0.8001891
    ##  [818,] 0.8080195
    ##  [819,] 1.0014790
    ##  [820,] 0.9593707
    ##  [821,] 0.8313592
    ##  [822,] 1.2167857
    ##  [823,] 1.0748758
    ##  [824,] 1.2167857
    ##  [825,] 1.0008314
    ##  [826,] 1.0252180
    ##  [827,] 1.1002390
    ##  [828,] 1.2768584
    ##  [829,] 1.3379525
    ##  [830,] 1.2167857
    ##  [831,] 1.0014790
    ##  [832,] 1.0748758
    ##  [833,] 0.9715200
    ##  [834,] 1.1002390
    ##  [835,] 0.8499509
    ##  [836,] 0.8080195
    ##  [837,] 1.1002390
    ##  [838,] 1.0008314
    ##  [839,] 0.8646186
    ##  [840,] 1.1372653
    ##  [841,] 0.8386694
    ##  [842,] 0.8080195
    ##  [843,] 0.8001891
    ##  [844,] 1.1372653
    ##  [845,] 0.8804854
    ##  [846,] 1.1372653
    ##  [847,] 1.0008314
    ##  [848,] 0.9715200
    ##  [849,] 1.5518528
    ##  [850,] 1.3379525
    ##  [851,] 1.1641895
    ##  [852,] 0.8080195
    ##  [853,] 1.0735319
    ##  [854,] 1.2768584
    ##  [855,] 1.1002390
    ##  [856,] 1.1233072
    ##  [857,] 1.2768584
    ##  [858,] 1.2167857
    ##  [859,] 1.0440420
    ##  [860,] 1.1233072
    ##  [861,] 0.9349355
    ##  [862,] 0.8313592
    ##  [863,] 0.9107817
    ##  [864,] 1.1002390
    ##  [865,] 0.8646186
    ##  [866,] 1.0252180
    ##  [867,] 1.3379525
    ##  [868,] 0.8759996
    ##  [869,] 0.8313592
    ##  [870,] 1.0735319
    ##  [871,] 0.9715200
    ##  [872,] 1.5518528
    ##  [873,] 1.0498496
    ##  [874,] 1.1002390
    ##  [875,] 0.8001891
    ##  [876,] 1.1233072
    ##  [877,] 0.9107817
    ##  [878,] 1.0498496
    ##  [879,] 1.0440420
    ##  [880,] 0.8804854
    ##  [881,] 1.0440420
    ##  [882,] 1.1233072
    ##  [883,] 0.8901779
    ##  [884,] 0.9390040
    ##  [885,] 1.2768584
    ##  [886,] 0.8759996
    ##  [887,] 0.8856563
    ##  [888,] 1.0748758
    ##  [889,] 1.0252180
    ##  [890,] 0.8856563
    ##  [891,] 0.8499509
    ##  [892,] 1.1486690
    ##  [893,] 0.8856563
    ##  [894,] 0.9121887
    ##  [895,] 1.1233072
    ##  [896,] 1.0498496
    ##  [897,] 1.0440420
    ##  [898,] 1.3379525
    ##  [899,] 1.0008314
    ##  [900,] 0.9349355
    ##  [901,] 1.0498496
    ##  [902,] 0.9715200
    ##  [903,] 1.2768584
    ##  [904,] 0.9715200
    ##  [905,] 0.8901779
    ##  [906,] 1.1641895
    ##  [907,] 0.8313592
    ##  [908,] 0.9349355
    ##  [909,] 0.8949477
    ##  [910,] 0.8386694
    ##  [911,] 0.8001891
    ##  [912,] 0.8804854
    ##  [913,] 1.1641895
    ##  [914,] 1.5518528
    ##  [915,] 0.9349355
    ##  [916,] 0.8646186
    ##  [917,] 1.1641895
    ##  [918,] 1.0440420
    ##  [919,] 0.8499509
    ##  [920,] 1.0748758
    ##  [921,] 0.9390040
    ##  [922,] 0.9349355
    ##  [923,] 1.0440420
    ##  [924,] 1.1641895
    ##  [925,] 0.8646186
    ##  [926,] 1.1002390
    ##  [927,] 0.8856563
    ##  [928,] 0.8646186
    ##  [929,] 1.1486690
    ##  [930,] 0.9121887
    ##  [931,] 1.0440420
    ##  [932,] 1.1233072
    ##  [933,] 0.9593707
    ##  [934,] 0.9107817
    ##  [935,] 1.1002390
    ##  [936,] 1.1233072
    ##  [937,] 1.0748758
    ##  [938,] 0.8499509
    ##  [939,] 1.2768584
    ##  [940,] 0.8646186
    ##  [941,] 0.8901779
    ##  [942,] 0.8499509
    ##  [943,] 1.1641895
    ##  [944,] 1.2768584
    ##  [945,] 1.0014790
    ##  [946,] 1.2167857
    ##  [947,] 0.8646186
    ##  [948,] 0.9593707
    ##  [949,] 0.8856563
    ##  [950,] 0.9390040
    ##  [951,] 0.8595390
    ##  [952,] 1.5518528
    ##  [953,] 1.1372653
    ##  [954,] 0.8856563
    ##  [955,] 0.8386694
    ##  [956,] 0.8499509
    ##  [957,] 1.0014790
    ##  [958,] 1.0498496
    ##  [959,] 1.1486690
    ##  [960,] 0.8499509
    ##  [961,] 0.9349355
    ##  [962,] 0.8386694
    ##  [963,] 1.1002390
    ##  [964,] 0.8386694
    ##  [965,] 1.1486690
    ##  [966,] 0.9715200
    ##  [967,] 1.2768584
    ##  [968,] 1.1233072
    ##  [969,] 0.9593707
    ##  [970,] 0.9349355
    ##  [971,] 1.0498496
    ##  [972,] 0.8968552
    ##  [973,] 0.9390040
    ##  [974,] 1.3379525
    ##  [975,] 0.8386694
    ##  [976,] 0.8313592
    ##  [977,] 0.9593707
    ##  [978,] 1.1002390
    ##  [979,] 0.9593707
    ##  [980,] 0.9390040
    ##  [981,] 0.8901779
    ##  [982,] 1.2768584
    ##  [983,] 0.9390040
    ##  [984,] 0.9349355
    ##  [985,] 1.0440420
    ##  [986,] 1.0252180
    ##  [987,] 1.1372653
    ##  [988,] 0.8386694
    ##  [989,] 0.8901779
    ##  [990,] 1.0008314
    ##  [991,] 0.8386694
    ##  [992,] 0.8595390
    ##  [993,] 1.2167857
    ##  [994,] 0.8968552
    ##  [995,] 0.9715200
    ##  [996,] 1.2768584
    ##  [997,] 0.8001891
    ##  [998,] 0.9593707
    ##  [999,] 1.0014790
    ## 
    ## $model.matrix
    ##   (Intercept) site_pH1
    ## 1           1        1
    ## 2           1        1
    ## 3           1        1
    ## 4           1        1
    ## 5           1       -1
    ## 6           1       -1
    ## 7           1       -1
    ## 8           1       -1
    ## 
    ## $terms
    ## dist_tab_assay ~ site_pH
    ## attr(,"variables")
    ## list(dist_tab_assay, site_pH)
    ## attr(,"factors")
    ##                site_pH
    ## dist_tab_assay       0
    ## site_pH              1
    ## attr(,"term.labels")
    ## [1] "site_pH"
    ## attr(,"order")
    ## [1] 1
    ## attr(,"intercept")
    ## [1] 1
    ## attr(,"response")
    ## [1] 1
    ## attr(,".Environment")
    ## <environment: R_GlobalEnv>
    ## 
    ## attr(,"class")
    ## [1] "adonis"

``` r
anova(betadisper(dist_tab_assay,samplesNatSim$site_pH))
```

    ## Analysis of Variance Table
    ## 
    ## Response: Distances
    ##           Df Sum Sq Mean Sq F value Pr(>F)
    ## Groups     1 306.92 306.921  3.5658 0.1079
    ## Residuals  6 516.45  86.074

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
    ##  [1] limma_3.50.3                EnhancedVolcano_1.12.0     
    ##  [3] ashr_2.2-54                 apeglm_1.16.0              
    ##  [5] tximport_1.22.0             ggvenn_0.1.9               
    ##  [7] dplyr_1.0.9                 vegan_2.6-2                
    ##  [9] lattice_0.20-45             permute_0.9-7              
    ## [11] gplots_3.1.3                genefilter_1.76.0          
    ## [13] RColorBrewer_1.1-3          markdown_1.1               
    ## [15] ggrepel_0.9.1               ggplot2_3.3.6              
    ## [17] BiocManager_1.30.18         devtools_2.4.3             
    ## [19] usethis_2.1.6               DESeq2_1.34.0              
    ## [21] SummarizedExperiment_1.24.0 Biobase_2.54.0             
    ## [23] MatrixGenerics_1.6.0        matrixStats_0.62.0         
    ## [25] GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
    ## [27] IRanges_2.28.0              S4Vectors_0.32.4           
    ## [29] BiocGenerics_0.40.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] plyr_1.8.7             splines_4.1.2          BiocParallel_1.28.3   
    ##   [4] digest_0.6.29          invgamma_1.1           htmltools_0.5.2       
    ##   [7] SQUAREM_2021.1         fansi_1.0.3            magrittr_2.0.3        
    ##  [10] memoise_2.0.1          cluster_2.1.3          tzdb_0.3.0            
    ##  [13] remotes_2.4.2          Biostrings_2.62.0      readr_2.1.2           
    ##  [16] annotate_1.72.0        extrafont_0.18         vroom_1.5.7           
    ##  [19] extrafontdb_1.0        bdsmatrix_1.3-6        prettyunits_1.1.1     
    ##  [22] colorspace_2.0-3       blob_1.2.3             xfun_0.31             
    ##  [25] callr_3.7.0            crayon_1.5.1           RCurl_1.98-1.6        
    ##  [28] survival_3.3-1         glue_1.6.2             gtable_0.3.0          
    ##  [31] zlibbioc_1.40.0        XVector_0.34.0         DelayedArray_0.20.0   
    ##  [34] proj4_1.0-11           pkgbuild_1.3.1         Rttf2pt1_1.3.10       
    ##  [37] maps_3.4.0             scales_1.2.0           mvtnorm_1.1-3         
    ##  [40] DBI_1.1.2              Rcpp_1.0.8.3           xtable_1.8-4          
    ##  [43] emdbook_1.3.12         bit_4.0.4              truncnorm_1.0-8       
    ##  [46] httr_1.4.3             ellipsis_0.3.2         farver_2.1.0          
    ##  [49] pkgconfig_2.0.3        XML_3.99-0.9           locfit_1.5-9.5        
    ##  [52] utf8_1.2.2             labeling_0.4.2         tidyselect_1.1.2      
    ##  [55] rlang_1.0.2            AnnotationDbi_1.56.2   munsell_0.5.0         
    ##  [58] tools_4.1.2            cachem_1.0.6           cli_3.3.0             
    ##  [61] generics_0.1.2         RSQLite_2.2.14         evaluate_0.15         
    ##  [64] stringr_1.4.0          fastmap_1.1.0          yaml_2.3.5            
    ##  [67] processx_3.5.3         knitr_1.39             bit64_4.0.5           
    ##  [70] fs_1.5.2               caTools_1.18.2         purrr_0.3.4           
    ##  [73] KEGGREST_1.34.0        nlme_3.1-157           ash_1.0-15            
    ##  [76] ggrastr_1.0.1          brio_1.1.3             compiler_4.1.2        
    ##  [79] rstudioapi_0.13        beeswarm_0.4.0         curl_4.3.2            
    ##  [82] png_0.1-7              testthat_3.1.4         tibble_3.1.7          
    ##  [85] geneplotter_1.72.0     stringi_1.7.6          ps_1.7.0              
    ##  [88] desc_1.4.1             ggalt_0.4.0            Matrix_1.4-1          
    ##  [91] vctrs_0.4.1            pillar_1.7.0           lifecycle_1.0.1       
    ##  [94] bitops_1.0-7           irlba_2.3.5            R6_2.5.1              
    ##  [97] KernSmooth_2.23-20     vipor_0.4.5            sessioninfo_1.2.2     
    ## [100] MASS_7.3-57            gtools_3.9.2.1         assertthat_0.2.1      
    ## [103] pkgload_1.2.4          rprojroot_2.0.3        withr_2.5.0           
    ## [106] GenomeInfoDbData_1.2.7 mgcv_1.8-40            parallel_4.1.2        
    ## [109] hms_1.1.1              coda_0.19-4            rmarkdown_2.14        
    ## [112] mixsqp_0.3-43          bbmle_1.0.25           numDeriv_2016.8-1.1   
    ## [115] ggbeeswarm_0.6.0
