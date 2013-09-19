#/bin/bash
##
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Copyright (C) 2013 Memorial Sloan-Kettering Cancer Center 
##

set -e 

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/dexseq_config.sh

echo
echo ${PROG}: Oqtans http://galaxy.cbio.mskcc.org Galaxy wrapper for the DEXSeq version ${DEXSEQ_VERSION}.
echo
echo DEXSeq: Detecting differential usage of exons from RNA-seq data.
echo 

## input arguments from the interface 
GFF_IN=${1}
shift
MATE_PAIR=${1}
shift
LIBTP=${1}
shift
minQL=${1}
shift
RES_FILE=${1}
shift
RES_WD=${1}
shift

## associated array with sequencing type.
declare -A SEQ_TYPE=( [no]=SE [yes]=PE )

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo

mkdir -p ${RES_WD}
echo extra file path $RES_WD
tmpGTF=`mktemp --tmpdir=/tmp`

echo load the genome annotation in GFF file 

${PYTHON_PATH} ${DIR}/dexseq_prepare_annotation.py ${GFF_IN} ${tmpGTF}
echo genome annotation stored in ${tmpGTF} 
echo

echo %%%%%%%%%%%%%%%%%%%%
echo % 2. Read counting %
echo %%%%%%%%%%%%%%%%%%%%
echo

tmpFILE=`mktemp --tmpdir=/tmp`
echo $tmpFILE
echo -e '\t'condition'\t'libType > ${tmpFILE}_CONDITIONS.tab

COND=0
for REPLICATE_GROUP in $@
do
    IFS=':'
    COND=$((COND+1))
    for BAM_FILE in ${REPLICATE_GROUP}
    do
        ## different group information 
        REPNAME=$(basename ${BAM_FILE%.dat})
        echo -e ${REPNAME}"\t"$COND"\t"${SEQ_TYPE[$MATE_PAIR]} >> ${tmpFILE}_CONDITIONS.tab
        
        ## counting the reads 
        ${SAMTOOLS_DIR}/samtools view -h $BAM_FILE | ${PYTHON_PATH} ${DIR}/dexseq_count.py -p ${MATE_PAIR} -s ${LIBTP} -a ${minQL} ${tmpGTF} - ${RES_WD}/${REPNAME}

        echo 
    done
    echo conuted condition ${COND} 
done
echo counted reads map to each exon.
echo 

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo "cat ${DIR}/run_DEXseq.R | $R_PATH --slave --args $tmpFILE $RES_WD $tmpGTF ${RES_FILE} $#" 
cat ${DIR}/run_DEXseq.R | $R_PATH --slave --args $tmpFILE $RES_WD $tmpGTF ${RES_FILE} 

## clean up
rm -fr ${RES_WD} ${tmpGTF} ${tmpFILE} 
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
