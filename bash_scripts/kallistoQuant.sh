#!/bin/bash

<<Block_comment 
kallistoQuant : Run the Kallisto tool on the data specified as arguments when the script is launched. 
These arguments are the type of organism (juvenile or adult), the type of reference transcriptome (A, LA, LJA), 
the type of sequences (paired or single), the name of dataset directory. Also, if the sequence type is single end, 
the average length of the RNA fragments (length) as well as the standard deviation of this metric (SD) must be specified. 
They can be easily obtained using the Python script averageLength_SD_bioanalyser.

Marc Meynadier
Block_comment

ORG_TYPE=$1 # adult or juvenile
TRANSCRIPTOME_TYPE=$2 # A (adult) or LA (Larvae Adult) or LJA (Larvae Juvenile Adult)
SEQUENCE_TYPE=$3 # paired or single (paired end or single end)

if [ $ORG_TYPE = "adult" ]; then
	SEQ_PATH=/scratch2/genomes/mmeynadier/results/trimmomatic/$ORG_TYPE
	SEQ_DIR=$4 # name of directory
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
