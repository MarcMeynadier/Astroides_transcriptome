#!/bin/bash

ORG_TYPE=$1 # adult or juvenile
TRANSCRIPTOME_TYPE=$2 # A (adult) or LA (Larvae Adult) or LJA (Larvae Juvenile Adult)
SEQUENCE_TYPE=$3 # paired or single

if [ $ORG_TYPE = "adult" ]; then
	SEQ_PATH=/scratch2/genomes/mmeynadier/results/trimmomatic/$ORG_TYPE
	SEQ_DIR=$4 # name of dir
fi

if [ $ORG_TYPE = "juvenile" ]; then
        SEQ_PATH=/backup2/genomes/mmeynadier/rawSequences/juvenile
fi

if [ $TRANSCRIPTOME_TYPE = "A" ]; then
  TRANSCRIPT_INDEX=/backup2/genomes/mmeynadier/kallisto/transcriptAdultIndex.idx
  OUTPUT_DIR=/scratch2/genomes/mmeynadier/results/kallisto/$ORG_TYPE/adultTranscriptome/$SEQ_DIR
fi

if [ $TRANSCRIPTOME_TYPE = "LA" ] && [ $ORG_TYPE = "adult" ]; then
  TRANSCRIPT_INDEX=/backup2/genomes/mmeynadier/kallisto/transcriptLarvaeAdultIndex.idx
  OUTPUT_DIR=/scratch2/genomes/mmeynadier/results/kallisto/$ORG_TYPE/larvaeAdultTranscriptome/$SEQ_DIR
fi

if [ $TRANSCRIPTOME_TYPE = "LA" ] && [ $ORG_TYPE = "juvenile" ]; then
  TRANSCRIPT_INDEX=/backup2/genomes/mmeynadier/kallisto/transcriptLarvaeAdultIndex.idx
  OUTPUT_DIR=/scratch2/genomes/mmeynadier/results/kallisto/$ORG_TYPE/larvaeAdultTranscriptome
fi

if [ $TRANSCRIPTOME_TYPE = "LJA" ] && [ $ORG_TYPE = "adult" ]; then
  TRANSCRIPT_INDEX=/backup2/genomes/mmeynadier/kallisto/transcriptLarvaeJuvenileAdultIndex.idx
  OUTPUT_DIR=/backup2/genomes/mmeynadier/kallisto//larvaeJuvenileAdultTranscriptome/$ORG_TYPE/$SEQ_DIR
fi

if [ $TRANSCRIPTOME_TYPE = "LJA" ] && [ $ORG_TYPE = "juvenile" ]; then
  TRANSCRIPT_INDEX=/backup2/genomes/mmeynadier/kallisto/transcriptLarvaeJuvenileAdultIndex.idx
  OUTPUT_DIR=/backup2/genomes/mmeynadier/kallisto/larvaeJuvenileAdultTranscriptome/$ORG_TYPE
fi

if [ $SEQUENCE_TYPE = "paired" ] && [ $ORG_TYPE = "adult" ]; then
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

if [ $SEQUENCE_TYPE = "paired" ] && [ $ORG_TYPE = "juvenile" ]; then
  for i in ${SEQ_PATH}/*_1.fq.gz
  do
    SPL_NAME=${i//_1/_}
    SAMPLE_NAME=$(echo $SPL_NAME | cut -d . -f 1 | sed 's|.*/||')
    j=${i//_1/_2}
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
