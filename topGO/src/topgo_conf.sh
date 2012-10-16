#!/bin/bash
## PATH set-up
#export PYTHOPATH

## get session path and result file 
E_FILE_PATH=${1}
shift
RES_FILE=${1}
shift 

## run topGO under extra file path and move the result file
mkdir -p $E_FILE_PATH
cd ${E_FILE_PATH}
python ${OQTANS_PATH}/topGO/src/topgo_main.py $@ > logfile.txt 
mv resultFile* $RES_FILE

## clean up 
rm -fr ${E_FILE_PATH}
## done 
