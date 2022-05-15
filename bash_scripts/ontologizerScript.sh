#!/bin/bash

<<Block_comment 
Ontologizer Script : Runs Ontologizer tool 

Marc Meynadier
Block_comment


pushd ../..//data/net/8_functionalAnnotation/ontologizer
java -jar Ontologizer.jar -a associationFile.ids -g go.obo -s studySamples -p populationFile.txt -c Parent-Child-Union -o outputResults/codingDEG -d 0.05 -r 1000
popd