#!/bin/bash

SEQ_PATH=$1
SEQ_PATH2=$2
SEQ_DIR=$3
OUTPUT_DIR=/scratch2/genomes/mmeynadier/results/md5Hashes

md5deep -rel $SEQ_PATH/* > $OUTPUT_DIR/InputHashes${SEQ_DIR}.md5
md5deep -rel $SEQ_PATH2/* > $OUTPUT_DIR/OutputHashes${SEQ_DIR}.md5

sort $OUTPUT_DIR/InputHashes${SEQ_DIR}.md5 > $OUTPUT_DIR/sortedInputHashes${SEQ_DIR}.txt
sort $OUTPUT_DIR/OutputHashes${SEQ_DIR}.md5 > $OUTPUT_DIR/sortedOutputHashes${SEQ_DIR}.txt

md5deep $OUTPUT_DIR/sortedInputHashes${SEQ_DIR}.txt > $OUTPUT_DIR/HashesComparisons${SEQ_DIR}.txt
md5deep $OUTPUT_DIR/sortedOutputHashes${SEQ_DIR}.txt >> $OUTPUT_DIR/HashesComparisons${SEQ_DIR}.txt

cat $OUTPUT_DIR/HashesComparisons${SEQ_DIR}.txt

