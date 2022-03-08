#!/bin/bash

ORG_TYPE=$1 # adult or juvenile
TRANSCRIPTOME_TYPE=$2 # A (adult) or LA (Larvae Adult) or LJA (Larvae Juvenile Adult)
SEQ_PATH=/scratch2/genomes/mmeynadier/results/trimmomatic/$ORG_TYPE
SEQ_DIR=$3 # name of dir
SEQUENCE_TYPE=$4 # paired or single

if [ $TRANSCRIPTOME_TYPE = "A" ]; then
  TRANSCRIPT_INDEX=/scratch2/genomes/mmeynadier/results/kallisto/transcriptAdultIndex.idx
  OUTPUT_DIR=/scratch2/genomes/mmeynadier/results/kallisto/$ORG_TYPE/adultTranscriptome/$SEQ_DIR
fi

if [ $TRANSCRIPTOME_TYPE = "LA" ]; then
  TRANSCRIPT_INDEX=/scratch2/genomes/mmeynadier/results/kallisto/transcriptLarvaeAdultIndex.idx
  OUTPUT_DIR=/scratch2/genomes/mmeynadier/results/kallisto/$ORG_TYPE/larvaeAdultTranscriptome/$SEQ_DIR
fi


if [ $SEQUENCE_TYPE = "paired" ]; then
  for i in ${SEQ_PATH}/${SEQ_DIR}/*_1_paired_trimmed.fq.gz
  do
    SPL_NAME=${i//_1_/_}
    SAMPLE_NAME=$(echo $SPL_NAME | cut -d . -f 1 | sed 's|.*/||')
    j=${i//_1_/_2_}
    kallisto quant -i $TRANSCRIPT_INDEX -o $OUTPUT_DIR -b 100 -t 10 $i $j
    mv ${OUTPUT_DIR}/abundance.h5 ${OUTPUT_DIR}/abundance_${SAMPLE_NAME}.h5 
    mv ${OUTPUT_DIR}/abundance.tsv ${OUTPUT_DIR}/abundance_${SAMPLE_NAME}.tsv
    mv ${OUTPUT_DIR}/run_info.json ${OUTPUT_DIR}/run_info_${SAMPLE_NAME}.json 
  done
fi

if [ $SEQUENCE_TYPE = "single" ]; then
  LENGTH=$5
  SD=$6
  for i in ${SEQ_PATH}/${SEQ_DIR}/*.fq.gz
  do
    SAMPLE_NAME=$(echo $i | cut -d . -f 1 | sed 's|.*/||')
    kallisto quant -i $TRANSCRIPT_INDEX -o $OUTPUT_DIR -b 100 -t 10 --single -l $LENGTH -s $SD $i
    mv ${OUTPUT_DIR}/abundance.h5 ${OUTPUT_DIR}/abundance_${SAMPLE_NAME}.h5
    mv ${OUTPUT_DIR}/abundance.tsv ${OUTPUT_DIR}/abundance_${SAMPLE_NAME}.tsv
    mv ${OUTPUT_DIR}/run_info.json ${OUTPUT_DIR}/run_info_${SAMPLE_NAME}.json
  done
fi
