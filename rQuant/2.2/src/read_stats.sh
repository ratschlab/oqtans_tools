#!/bin/bash

#
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
echo ReadStats generates a statistic about the read alignments
echo and the covered genes.
echo 

if [ -z "$11" -o "$1" == '--help' ];
then
  echo Usage: $0
  echo "   or:" $0 --help
  false
fi 

GFF3_INPUT=${1}
BAM_INPUT=${2}
RQUANT_RES_FILE=${3}
RQUANT_RES_DIR=${4}
GENES_FN=${RQUANT_RES_DIR}/genes.mat

mkdir -p $RQUANT_RES_DIR

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo
 
echo load the genome annotation in GFF3 format and create an annotation object
if [ ! -f ${GENES_FN} ]
then
    export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
    ${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${GFF3_INPUT} all all ${GENES_FN}
    ${DIR}/../bin/genes_cell2struct ${GENES_FN}
fi

echo
echo %%%%%%%%%%%%%%%%%%%%%
echo % 2. Read statistic %
echo %%%%%%%%%%%%%%%%%%%%%
echo

echo generate statistic about read alignments and covered genes
${DIR}/../bin/read_stats ${RQUANT_RES_DIR} ${BAM_INPUT} ${RQUANT_RES_FILE}

echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
