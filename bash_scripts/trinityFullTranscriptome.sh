#!/bin/bash

<<Block_comment 
trinityFullTranscriptome: Run the Trinity tool on the data specified as arguments in the code. 
The arguments are the type of sequences (fq, i.e. fastq), the maximum memory allocated, the output report in txt format, 
the computational power allocated, the output directory of the reference transcriptome as well as the argument full_cleanup, 
allowing to purge the directories of intermediate results created by the algorithm

Marc Meynadier
Block_comment

/scratch1/genomes/copley/src/trinityrnaseq-v2.13.2/Trinity --seqType fq --max_memory 128G --samples_file samples_file.txt --CPU 32 --output full_trinity --full_cleanup


