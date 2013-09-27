#!/bin/bash

set -e 

source ../../oqtans_conf.sh

cd $OQTANS_PATH/palmapper/0.5
if [ "$1" == "" -o "$1" == "all" ];
then
    make all
fi

if [ "$1" == "clean" ];
then
    make clean
fi
