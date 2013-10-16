#!/bin/bash
#
# oqtans module compile script
#
set -e 

source ./../../../oqtans_config.sh

cd $OQTANS_PATH/Trinity/r2013_08_14
if [ "$1" == "" -o "$1" == "all" ];
then
    make 
fi

if [ "$1" == "clean" ];
then
    make clean
fi
