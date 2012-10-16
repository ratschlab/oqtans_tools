#!/bin/bash

set -e 

source ../../oqtans_conf.sh

cd $OQTANS_PATH/DESeq/1.6/mex
if [ "$1" == "" -o "$1" == "all" ];
then
    make octave
fi

if [ "$1" == "clean" ];
then
    make clean
fi
