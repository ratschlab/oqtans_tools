#!/bin/bash
#
# oqtans submodule compile program
# 
set -e 

source ./../../../oqtans_config.sh

cd $OQTANS_PATH/BWA/0.6.2
if [ "$1" == "" -o "$1" == "all" ];
then
    make 
fi

if [ "$1" == "clean" ];
then
    make clean
fi
