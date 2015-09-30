#!/bin/bash
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
#

set -e 

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/rquant_config.sh

echo
echo ${PROG}: This program is part of the rQuant version $RQUANT_VERSION.
echo
echo rQuant determines the abundance of multiple transcripts per 
echo gene locus from RNA-Seq measurements.
echo 

ANNO_INPUT=${1}
BAM_INPUT=${2}
RQUANT_RES_FILE=${3}
RQUANT_RES_DIR=${4}
LOAD_PROFILES=${5}
PROFILES_FN=${6}
LEARN_PROFILES=${7}
PROFILES_FN_OUT=${8}

mkdir -p $RQUANT_RES_DIR
touch $RQUANT_RES_DIR/genes.mat

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo
echo load the genome annotation in GFF3 format and create an annotation object
echo
export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
${PYTHON_PATH} ${DIR}/../tools/GFFParser.py ${ANNO_INPUT} ${RQUANT_RES_DIR}/genes.mat

echo
echo %%%%%%%%%%%%%%%%%%%%%
echo % 2. Quantification %
echo %%%%%%%%%%%%%%%%%%%%%
echo

echo quantify transcripts using given alignments
echo 
(${DIR}/../bin/rquant ${RQUANT_RES_DIR} ${BAM_INPUT} ${RQUANT_RES_FILE} ${RQUANT_RES_DIR}/ ${LOAD_PROFILES} ${PROFILES_FN} ${LEARN_PROFILES} ${PROFILES_FN_OUT} 2>&1 || (echo rquant failed 1>&2))
echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
