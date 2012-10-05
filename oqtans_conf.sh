#!/bin/bash

export OQTANS_ROOT_PATH=/mnt/oqtansTools
export OQTANS_PATH=$OQTANS_ROOT_PATH/oqtans
export OQTANS_SRC_PATH=$OQTANS_ROOT_PATH/oqtans_src
export OQTANS_DEP_PATH=$OQTANS_ROOT_PATH/oqtans_dep

export OQTANS_TMP_PATH=$OQTANS_ROOT_PATH/tmp
mkdir -p $OQTANS_TMP_PATH

export OQTANS_GIT_BRANCH=ubuntu-10.04

if [ "$1" == "-s" ];
then
    echo $OQTANS_SRC_PATH
fi
