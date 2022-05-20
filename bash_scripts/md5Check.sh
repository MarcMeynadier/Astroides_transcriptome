#!/bin/bash

<<Block_comment 
md5Check : Run the md5deep tool which allows to obtain md5 hashes recursively on a whole directory. 
The hashes are then sorted and a hash per folder is obtained so that the user can compare this hash with another one. 
The arguments needed to launch the script are the path to the first folder to be compared, the path to the second folder to be compared, 
and the name of the dataset from which the compared data comes.

Marc Meynadier
Block_comment

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

