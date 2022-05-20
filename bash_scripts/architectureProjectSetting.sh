#!/bin/bash

<<Block_comment 
architectureProjectSettings : Construction and layout of the architecture of the specific folders for the Astroides_transcriptome project

Marc Meynadier
Block_comment

mkdir Astroides_transcriptome
cd Astroides_transcriptome

mkdir output data scripts

cd output
mkdir DESeq2 functionalGenesAnalysis

cd ../scripts
git clone https://github.com/MarcMeynadier/Astroides_transcriptome
mv Astroides_transcriptome/* .
rm -rf Astroides_transcriptome

cd ../data
mkdir raw net
cd raw
mkdir adult juvenile
mkdir -p adult/{rawSequences,kallisto}/{nov2016,june2017,sept2017,may2018,sept2018,specialSequences}

cd ../net 
mkdir 1_md5checksum 2_fastqcPreTrim 3_trimmomatic 4_fastqcPostTrim 5_trinity 6_kallisto 7_deseq2 8_functionalAnnotation
mkdir -p 2_fastqcPreTrim/{adult,juvenile}
mkdir -p 2_fastqcPreTrim/adult//{nov2016,june2017,sept2017,may2018,sept2018,specialSequences}
mkdir -p 3_trimmomatic/{adult,juvenile}
mkdir -p 3_trimmomatic/adult/{nov2016,june2017,sept2017,may2018,sept2018,specialSequences}
mkdir -p 4_fastqcPostTrim/{adult,juvenile}
mkdir -p 4_fastqcPostTrim/adult/{nov2016,june2017,sept2017,may2018,sept2018,specialSequences}
mkdir -p 6_kallisto/{adult,juvenile}
mkdir -p 6_kallisto/adult/{1_preliminarySamples,2_spatialComparison,3_temporalComparison,4_trueTransplant,5_gardenShort}
mkdir -p 7_deseq2/adulTranscriptome
mkdir -p 7_deseq2/adulTranscriptome/{adult,juvenile}
mkdir -p 8_functionalAnnotation/ontologizer
mkdir -p 8_functionalAnnotation/functionalGenesAnalysis
mkdir -p 8_functionalAnnotation/functionalGenesAnalysis/DESeq2_analysis
mkdir -p 8_functionalAnnotation/functionalGenesAnalysis/ontologizer_analysis
mkdir -p 8_functionalAnnotation/functionalGenesAnalysis/DESeq2_X_ontologizer


