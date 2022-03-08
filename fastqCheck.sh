#!/bin/bash

ORG_TYPE=$1 # adult or juvenile
SEQ_PATH_PRE=/backup2/genomes/mmeynadier/rawSequences/$ORG_TYPE
SEQ_PATH_POST=/scratch2/genomes/mmeynadier/results/trimmomatic/$ORG_TYPE
SEQ_DIR=$2 # name of dir
FORMAT=$3 # txt or fastq
TRIM_TYPE=$4 # preTrim or postTrim
OUTPUT_DIR_PRE=/scratch2/genomes/mmeynadier/results/fastqc_preTrimming/$ORG_TYPE
OUTPUT_DIR_POST=/scratch2/genomes/mmeynadier/results/fastqc_postTrimming/$ORG_TYPE

if [ $TRIM_TYPE = "preTrim" ]; then
  fastqc -o $OUTPUT_DIR_PRE/$SEQ_DIR $SEQ_PATH_PRE/$SEQ_DIR/*${FORMAT}
fi

if [ $TRIM_TYPE = "postTrim" ]; then
  fastqc -o $OUTPUT_DIR_POST/$SEQ_DIR $SEQ_PATH_POST/$SEQ_DIR/*${FORMAT}
fi

