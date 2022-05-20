#!/bin/bash

<<Block_comment 
fastqCheck : Run the FastQC tool on the data specified as arguments when the script is launched. 
These arguments are the type of organism (juvenile or adult), the format of the sequences (txt or fastq), 
as well as whether the sequences have already been trimmed or not (preTrim or postTrim). 
The results of FastQC are stored in a folder related to the results of this tool, specified in the project architecture.

Marc Meynadier
Block_comment

ORG_TYPE=$1 # adult or juvenile
SEQ_PATH_PRE=/backup2/genomes/mmeynadier/rawSequences/$ORG_TYPE
SEQ_PATH_POST=/scratch2/genomes/mmeynadier/results/trimmomatic/$ORG_TYPE
FORMAT=$2 # txt or fastq
TRIM_TYPE=$3 # preTrim or postTrim
OUTPUT_DIR_PRE=/scratch2/genomes/mmeynadier/results/fastqc_preTrimming/$ORG_TYPE
OUTPUT_DIR_POST=/scratch2/genomes/mmeynadier/results/fastqc_postTrimming/$ORG_TYPE

if [ $TRIM_TYPE = "preTrim" ]; then
  fastqc -o $OUTPUT_DIR_PRE $SEQ_PATH_PRE/*${FORMAT}
fi

if [ $TRIM_TYPE = "postTrim" ]; then
  fastqc -o $OUTPUT_DIR_POST $SEQ_PATH_POST/*${FORMAT}
fi

