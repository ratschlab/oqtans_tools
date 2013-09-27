#!/bin/bash
#
# oqtans sub-module compile script  
#
set -e 

if [ "$1" == "" -o "$1" == "all" ];
then
    make 
    cp mmr $OQTANS_DEP_PATH/bin/
fi
