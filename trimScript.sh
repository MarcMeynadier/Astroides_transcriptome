#!/bin/bash

ORG_TYPE=$1 # adult or juvenile
SEQ_PATH=/backup2/genomes/mmeynadier/rawSequences/$ORG_TYPE
SEQ_DIR=$2 # name of dir
FORMAT=$3 # txt or fastq
SEQUENCE_TYPE=$4 # paired or single
OUTPUT_DIR=/scratch2/genomes/mmeynadier/results/trimmomatic/$ORG_TYPE/$SEQ_DIR
TRIMMOMATIC=/usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar


if [ $SEQUENCE_TYPE = "paired" ]; then
  ADAPTERS=ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10
  for i in ${SEQ_PATH}/${SEQ_DIR}/*1.${FORMAT}
  do
    j=${i%%1.${FORMAT}}"2.${FORMAT}"
    SAMPLE_NAME=$(echo $i | cut -d . -f 1 | sed 's/\(.*\)_.*/\1/')
    java -jar $TRIMMOMATIC PE -threads 10 -phred33 -trimlog ${SAMPLE_NAME}_trimlog.txt \
    ${i} ${j} ${SAMPLE_NAME}_1_paired_trimmed.fq.gz ${SAMPLE_NAME}_1_unpaired_trimmed.fq.gz ${SAMPLE_NAME}_2_paired_trimmed.fq.gz ${SAMPLE_NAME}_2_unpaired_trimmed.fq.gz \
    $ADAPTERS LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  done

  mv ${SEQ_PATH}/${SEQ_DIR}/*trimmed.fq.gz ${OUTPUT_DIR}/
  mv ${SEQ_PATH}/${SEQ_DIR}/*trimlog.txt ${OUTPUT_DIR}/
fi


if [ $SEQUENCE_TYPE = "single" ]; then 
  ADAPTERS=ILLUMINACLIP:/usr/local/bin/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10
  for i in ${SEQ_PATH}/${SEQ_DIR}/*.${FORMAT}
  do
    SAMPLE_NAME=$(echo $i | cut -d . -f 1)
    java -jar $TRIMMOMATIC SE -threads 10 -phred33 -trimlog ${SAMPLE_NAME}_trimlog.txt \
    ${i} ${SAMPLE_NAME}_trimmed.fq.gz \
    $ADAPTERS LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  done

  mv ${SEQ_PATH}/${SEQ_DIR}/*trimmed.fq.gz ${OUTPUT_DIR}/
  mv ${SEQ_PATH}/${SEQ_DIR}/*trimlog.txt ${OUTPUT_DIR}/
fi

