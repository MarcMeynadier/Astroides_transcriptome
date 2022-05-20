#!/bin/bash

<<Block_comment 
ontologizerScript : The majority of these arguments are files that could be obtained directly and placed in the directory 
containing the Ontologizer tool on the command line thanks to the Python script enrichmentAnalysisParser, integrated into 
the functionalGenesAnalysis program. These arguments are an association file (link between the genes and their Gene Ontology codes), 
the Gene Ontology database file, a study samples file containing the differentially expressed genes found in each contrast performed by DESeq2, 
a population file containing all the genes contained in the reference transcriptome provided by Trinity, the calculation method of the algorithm 
as well as the output directory of the results

Marc Meynadier
Block_comment

pushd ../..//data/net/8_functionalAnnotation/ontologizer
java -jar Ontologizer.jar -a associationFile.ids -g go.obo -s studySamples -p populationFile.txt -c Parent-Child-Union -o outputResults/codingDEG -d 0.05 -r 1000
popd