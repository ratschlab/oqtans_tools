#!/bin/bash

set -e 

source ../../oqtans_conf.sh

cd $OQTANS_PATH/mTIM/0.2/src/utils
if [ "$1" == "" -o "$1" == "all" ];
then
    make octave
fi

if [ "$1" == "clean" ];
then
    make clean
fi
