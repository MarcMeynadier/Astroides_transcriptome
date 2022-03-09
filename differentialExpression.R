# Differential expression on Kallisto data 

setwd('~/Documents/Projet/code/kallistoResults')

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

packageCheckClassic(c('DESeq2','tidyverse','devtools','BiocManager','tximportData','rhdf5'))
#remotes::install_github("pachterlab/sleuth#260")
#BiocManager::install('tximport', force = TRUE)
#library('tximport')

# Data importation

BiocManager::install("tximportData")
library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)


adultTranscriptome_adult_nov2016 <- '~/Documents/Projet/code/kallistoResults/adultTranscriptome/adult/nov2016'
larvaeAdultTranscriptome_adult_nov2016 <- '~/Documents/Projet/code/kallistoResults/larvaeAdultdultTranscriptome/adult/nov2016'



