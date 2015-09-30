#!/bin/bash
#
# oqtans submodule compile program 
# 
set -e 

source ./../../../oqtans_config.sh

cd $OQTANS_PATH/rQuant/2.2/mex
if [ "$1" == "" -o "$1" == "all" ];
then
    make octave
fi

if [ "$1" == "clean" ];
then
    make clean
fi
