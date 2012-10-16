#!/bin/bash

set -e 

source ../../oqtans_conf.sh

cd $OQTANS_PATH/bwa/0.6.2
if [ "$1" == "" -o "$1" == "all" ];
then
    make 
fi

if [ "$1" == "clean" ];
then
    make clean
fi
