#/bin/bash
##
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Copyright (C) 2009-2013 Max Planck Society & Memorial Sloan-Kettering Cancer Center 
##

set -e 

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/edger_config.sh

echo
echo ${PROG}: Oqtans http://galaxy.cbio.mskcc.org Galaxy wrapper for the edgeR version $edgeR_VERSION.
echo
echo edgeR provides statistical routines for assessing differential expression in RNA-Seq experiments.
echo 

ANNO_GFF=${1}
shift
RES_FILE=${1}
shift
GENES_FN=${1}
shift

mkdir -p `dirname $GENES_FN`
touch ${GENES_FN}

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo
echo load the genome annotation and creating annotation object...
${PYTHON_PATH} ${DIR}/../tools/GFFParser.py ${ANNO_GFF} ${GENES_FN}
echo genome annotation object stored in $GENES_FN
echo 

ANL_TYPE=${1}
shift
FDR_MT=${1}
shift

echo %%%%%%%%%%%%%%%%%%%%
echo % 2. Read counting %
echo %%%%%%%%%%%%%%%%%%%%
echo
echo checking for BAM index... 
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
echo
echo counting reads overlapping exons using given alignments...
tmpfile=`mktemp --tmpdir=/tmp`

echo "${DIR}/../bin/get_read_counts ${GENES_FN} $tmpfile $@"
${DIR}/../bin/get_read_counts ${GENES_FN} $tmpfile "$@" 
echo counted reads map to each feature.
echo 

echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo
echo testing genes for differential expression using given read alignments

echo "cat ${DIR}/../src/edger_difftest.R | $R_PATH --slave --args $ANL_TYPE $FDR_MT $tmpfile ${RES_FILE}"
cat ${DIR}/../src/edger_difftest.R | $R_PATH --slave --args $ANL_TYPE $FDR_MT $tmpfile ${RES_FILE}  

echo %%%%%%%%
echo % Done %
echo %%%%%%%%
