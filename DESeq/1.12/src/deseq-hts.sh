#/bin/bash
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Copyright (C) 2009-2013 Max Planck Society & Memorial Sloan-Kettering Cancer Center 
#

set -e 

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/deseq_config.sh

echo
echo ${PROG}: Oqtans http://galaxy.cbio.mskcc.org Galaxy wrapper for the DESeq version $DESEQ_VERSION.
echo
echo DESeq performs differential expression testing from RNA-Seq measurements.
echo 

ANNO_INPUT=${1}
shift
DESEQ_RES_FILE=${1}
shift
GENES_FN=${1}
shift

mkdir -p `dirname $GENES_FN`

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo
echo load the genome annotation in GFF3 format and create an annotation object
export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${ANNO_INPUT} ${GENES_FN}
${DIR}/../bin/genes_cell2struct ${GENES_FN} 
echo 
echo genome annotation stored in $GENES_FN

echo %%%%%%%%%%%%%%%%%%%%
echo % 2. Read counting %
echo %%%%%%%%%%%%%%%%%%%%
echo

echo counting reads overlapping exons using given alignments
for REPLICATE_GROUP in $@
do
    IFS=':'
    for BAM_FILE in ${REPLICATE_GROUP}
    do
        echo
        if [ ! -f ${BAM_FILE}.bai ]
        then
            echo "Indexing $BAM_FILE"
            ${SAMTOOLS_DIR}/samtools index $BAM_FILE
        else
            echo "$BAM_FILE already indexed"
        fi
        echo
    done
done
tmpfile=`mktemp --tmpdir=/tmp`

echo "${DIR}/../bin/get_read_counts ${GENES_FN} $tmpfile $@"
${DIR}/../bin/get_read_counts ${GENES_FN} $tmpfile "$@" 

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo
echo testing genes for differential expression using given alignments

echo "cat ${DIR}/../src/difftest_deseq.R | $R_PATH --slave --args $tmpfile ${DESEQ_RES_FILE} $#"
cat ${DIR}/../src/difftest_deseq.R | $R_PATH --slave --args $tmpfile ${DESEQ_RES_FILE} $# 

rm $tmpfile ${tmpfile}_COUNTS.tab ${tmpfile}_CONDITIONS.tab
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
