#!/bin/bash
#
#
# 
set -e 

source ./../../../oqtans_config.sh

cd $OQTANS_PATH/EasySVM/0.3.3
if [ "$1" == "" -o "$1" == "all" ];
then
    python setup.py install --prefix=$OQTANS_DEP_PATH
fi

if [ "$1" == "clean" ];
then
    rm -rf build
fi
