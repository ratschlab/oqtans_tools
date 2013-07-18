#!/bin/bash

set -e 

source ../../oqtans_conf.sh

cd $OQTANS_PATH/DESeq/2-1.0.18/mex
if [ "$1" == "" -o "$1" == "all" ];
then
    make octave
fi

if [ "$1" == "clean" ];
then
    make clean
fi
