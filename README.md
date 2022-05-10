# Astroides_transcriptome

This repository is about comparative transcriptomic on the scleractinian coral Astroides calycularis.

Three different directories are different based on the script language, but these also correspond to the different stages of the comparative coral transcriptome study.

The first directory contains scripts written in bash, most of which are designed to be run on a remote server. The script architectureProjectSetting allows to create automatically the architecture (the set of folders and their layout) of the project, the majority of the scripts running locally being dependent on this architecture. md5Check allows to check the md5 hash of the raw data containing the sequences, in order to verify that these last ones are integrated before starting the pre-processing. fastqCheck allows to run the fastQC tool on raw sequences, either on raw sequences or on trim sequences after using the Trimmomatic tool. trimScript allows to run the Trimmomatic tool on raw sequences, taking into account the type of coral (juvenile or adult), the type of file (fasta or txt) and the type of sequencing (single end or paired end). trinityTranscriptome allows to launch the Trinity tool, taking into account an input text file containing the path to the sequences used to reconstruct the de novo transcriptome

