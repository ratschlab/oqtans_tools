#!/bin/bash

set -e 

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/rdiff_config.sh

echo
echo ${PROG}: This program is part of the rDiff version $RDIFF_VERSION.
echo
echo rDiff is a tool for differential expression testing of RNA-Seq data
echo

ANNO_INPUT=${1}
BAM_INPUT_1=${2}
BAM_INPUT_2=${3}
TEST_FN=${4}
RDIFF_RES_FILE=${5}
RDIFF_RES_DIR=${6}

mkdir -p $RDIFF_RES_DIR

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo
echo load the genome annotation in GFF3 format and create an annotation object
echo
export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${ANNO_INPUT} ${RDIFF_RES_DIR}/genes.mat
${DIR}/../bin/genes_cell2struct ${RDIFF_RES_DIR}/genes.mat 2>&1

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo % 2. Differential testing
echo %%%%%%%%%%%%%%%%%%%%%%%%%
echo
echo test method ${TEST_FN} 
echo 
if [ ! -f ${BAM_INPUT_1}.bai ]
then
	$SAMTOOLS_DIR/samtools index ${BAM_INPUT_1}
fi	
if [ ! -f ${BAM_INPUT_2}.bai ]
then
	$SAMTOOLS_DIR/samtools index ${BAM_INPUT_2}
fi	
${DIR}/../bin/galaxy_rdiff_web ${RDIFF_RES_DIR} ${BAM_INPUT_1} ${BAM_INPUT_2} ${RDIFF_RES_FILE} ${TEST_FN} 2>&1
echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
