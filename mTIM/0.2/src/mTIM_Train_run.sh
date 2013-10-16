#!/bin/bash
##
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Copyright (C) 2010-2013 Max Planck Society & Memorial Sloan-Kettering Cancer Center 
##

set -e 

DIR=`dirname $0`
PROG=`basename $0`

. ${DIR}/../bin/mTIM_config.sh

echo
echo ${PROG}: Oqtans http://galaxy.cbio.mskcc.org Galaxy wrapper for the mTIM version $MTIM_VERSION.
echo
echo mTIM margin-based Transcript Mapping from RNA-Sequencing data.
echo

GenomeSource=${1}
shift 
OUT_DIR=${1}
shift 

mkdir -p ${OUT_DIR}

echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '1. Genome Sequence Preparation'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo
if [ ${GenomeSource} = "fasta" ];
then
    FASTA_IN=${1}
    shift 
    echo Using the genome in fasta format ${FASTA_IN} 
    ${PYTHON_PATH} ${DIR}/../tools/GenomePrep.py $FASTA_IN $OUT_DIR 
    echo Genome information saved. 
else
    GIO_IN=${1}
    shift
fi 

GFF_IN=${1}
shift

BAM_IN=${1}
shift

ACC_SPF=${1}
shift
DON_SPF=${1}
shift

echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '2. Genome Annotation Preparation'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo

touch ${OUT_DIR}/gene.mat
${PYTHON_PATH} ${DIR}/../tools/GFFParser.py $GFF_IN $OUT_DIR/genes.mat

echo Genome annotation object saved.

echo '%%%%%%%%%%%'
echo '3. Training'
echo '%%%%%%%%%%%'
echo
if [ ! -f ${BAM_IN}.bai ]
then
	${SAMTOOLS_DIR}/samtools index ${BAM_IN}
fi

if [ ${GenomeSource} = "fasta" ];
then
    ${DIR}/../bin/mTIM_galaxy_train ${OUT_DIR} ${ACC_SPF} ${DON_SPF} ${BAM_IN} ${OUT_DIR} ${OUT_DIR} 
else
    ${DIR}/../bin/mTIM_galaxy_train ${GIO_IN} ${ACC_SPF} ${DON_SPF} ${BAM_IN} ${OUT_DIR} ${OUT_DIR} 
fi
echo '%%%%%%'
echo ' Done.'
echo '%%%%%%'
