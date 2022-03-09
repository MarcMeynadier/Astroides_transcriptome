# Differential expression on Kallisto data 

setwd('~/Documents/Projet/code/kallistoResults/adultTranscriptome/adult/nov2016')

# Functions

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

tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}

# Packages dependance

packageCheckClassic(c('DESeq2','tidyverse','devtools','BiocManager','tximport','rhdf5'))
#remotes::install_github("pachterlab/sleuth#260")
#BiocManager::install('tximport', force = TRUE)
#library('tximport')

# Data importation - tximport

library(tximport)

samples<-read.table('../sample_list.txt',header=T)

# regex try
grx <- glob2rx("*.tsv")

files <- file.path(samples$dir, with(samples, subset(samples, subset = grepl(grx, rownames(samples)))))

names(files)<-samples$sample

txi <- tximport(files, type = "kallisto", txOut=T)

#normal try
files <- file.path(samples$dir, "abundance.tsv")

names(files)<-samples$sample

txi <- tximport(files, type = "kallisto", txOut=T)

names(txi)

head(txi$counts)


adultTranscriptome_adult_nov2016 <- '~/Documents/Projet/code/kallistoResults/adultTranscriptome/adult/nov2016'
larvaeAdultTranscriptome_adult_nov2016 <- '~/Documents/Projet/code/kallistoResults/larvaeAdultdultTranscriptome/adult/nov2016'



