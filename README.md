# Astroides_transcriptome

<p align="center">
  <img src="https://github.com/MarcMeynadier/Astroides_transcriptome/blob/master/pictures/astroides.png?raw=true" alt="Astroides calycularis"/>
</p>
  
This repository is about comparative transcriptomic on the scleractinian coral *Astroides calycularis*.

Three different directories are different based on the script language, but these also correspond to the different stages of the comparative coral transcriptome study.

The first directory contains scripts written in bash, most of which are designed to be run on a remote server. 
* The script architectureProjectSetting allows to create automatically the architecture (the set of folders and their layout) of the project, the majority of the scripts running locally being dependent on this architecture. 
* md5Check allows to check the md5 hash of the raw data containing the sequences, in order to verify that these last ones are integrated before starting the pre-processing. 
* fastqCheck allows to run the fastQC tool on raw sequences, either on raw sequences or on trim sequences after using the Trimmomatic tool. 
* trimScript allows to run the Trimmomatic tool on raw sequences, taking into account the type of coral (juvenile or adult), the type of file (fasta or txt) and the type of sequencing (single end or paired end). 
* trinityTranscriptome allows to launch the Trinity tool, taking into account an input text file containing the path to the sequences used to reconstruct the de novo transcriptome. 
* Finally, kallistoQuant allows to run the Kallisto tool on both types of corals (adults and juveniles), three different types of transcriptomes having been built within this project, as well as the type of sequence (single end and paired end).

The second directory contains a set of R scripts specific to the use of the DESeq2 package, as well as a set of text files allowing the import of data from Kallisto using the tximport package. 
* The tx2gene files are also text files allowing to link a transcript to its gene of origin, and are issued from the reconstruction of the *de novo* transcriptome with the Trinity tool. Each R script contains the following structure: import of packages, data, pre-filtering of data, then calculation of the differential expression of genes. This differential expression is then expressed for a set of contrasts (comparison between two conditions) determined by the objectives of the study. After obtaining the differentially expressed genes for each contrast, they are visualized using MA-plot and volcano plot. Exploratory graphs are then made to explore the global results, *i.e.* a PCA as well as Venn diagrams, then a permanova is made to observe the potential statistical differences between the factors conditioning the data set. The results of the contrasts are then exported as CSV files.

The third directory is dedicated to Python scripts. These scripts are mostly dedicated to the exploitation of DESeq2 and Ontologizer results, except averageLength_SD_bioanalyzer. 
* averageLength_SD_bioanalyzer function is to obtain the average size of fragments (RNA or DNA) as well as the standard deviation of this size from the results of a Bioanalyzer instrument. Those informations is needed to run the bash script kallistoQuant in single end mode, and this script allows the user to get them quickly. 
* The other scripts belong to the functionalGenesAnalysis program, of which the eponymous script is the main script allowing the user to activate the other modules of the program. 
* The first module, DESeq2_analysis, allows the user to retrieve the results of the contrasts of differentially expressed genes, to filter them according to a p-value threshold, and to retrieve the functional annotations of the remaining genes. These annotations are of type Pfam, Gene Ontology as well as Panther. If the genes are coding, the protein sequences are also recovered and associated with the genes concerned. It is also possible to filter the results between several DESeq2 files, and to select the genes shared by these files, or conversely the genes that are not shared.
* The second module, Ontologizer_analysis, allows to retrieve the results from Ontologizer, to filter them according to a p-value threshold, and to format them as a dataframe. As with the previous module, it is possible to choose to sort the functionally enriched genes that are shared between different Ontologizer results, or those that are not shared. 
* The third module, DESeq2_X_ontologizer, is used to cross-reference the results obtained with the first two modules in order to retain only the genes in common between the two, for the same file type. If some genes are retained by the program, they are placed in a new CSV file. This module thus ensures that the user only retains differentially expressed genes whose Gene Ontology terms have been functionally enriched, *i.e.* whose presence is statistically significant. 
* The fourth module, finalExpressionAnalyzer, uses the results generated by the third module to associate one or more biological functions to each gene from its Gene Ontology terms. This sorting is performed thanks to a dictionary containing keywords specific to each biological function. From the output of this processing, a graph is generated thanks to the Matplotlib package and allows to display the cumulative expression value of the different functions from the log2 fold change of the different genes.
* The last module allows you to determine the settings of the program and to perform specific tasks, being connected to several scripts. The script enrichment_analyser_parser allows, from the files of the Pfam annotation (allowing to link the genes to their Gene ontology terms via the pfam2go database), the *de novo* transcriptome in output of Trinity as well as the results of DESeq2, to generate the various files necessary in input of Ontologizer. It is also possible to format the input files of the GO_MWU tool, another functional enrichment software, but which was not used afterwards. The threshold_settings script allows the user to choose the p-value used to filter the results of the first two modules, set by default to 0.05, as well as to activate or deactivate filtering by candidate genes. This last option is based on the last script, getCandidates, which generates the list of candidate genes through a protein keyword search on the Pfam annotation file of the reference transcriptome, which is less generic than the Gene Ontology terms. A CSV file is output containing the list of candidate genes. By activating the filtering by candidate genes, the program relies on this file to filter the results of the first two modules according to the genes concerned, ensuring robust results thanks to this particularly stringent method. 



