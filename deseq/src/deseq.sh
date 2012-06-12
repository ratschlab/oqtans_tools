#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch, Andre Kahles
# Copyright (C) 2009-2011 Max Planck Society
#

set -e 

PROG=`basename $0`
DIR=`dirname $0`

. ${DIR}/../bin/deseq_config.sh

echo
echo ${PROG}: This program is part of the DESeq version $DESEQ_VERSION.
echo
echo DESeq performs differential expression testing from RNA-Seq measurements.
echo 

if [ -z "$6" -o "$1" == '--help' ];
then
  echo Usage: $0
  echo "   or:" $0 --help
  false
fi 

ANNO_INPUT=${1}
shift
ANNO_FORMAT=${1}
shift
GENES_FN=${1}
shift
DESEQ_RES_FILE=${1}
shift


echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo

if [ "$ANNO_FORMAT" == '0' ]
then
    echo load the genome annotation in GFF3 format and create an annotation object
    echo
	IMfilename=${GENES_FN%.dat}.mat
    if [[ ! -f ${IMfilename} || ! -s ${IMfilename} ]]
    then
		export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
		${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${ANNO_INPUT} all all ${IMfilename}
    fi   
    ${DIR}/../bin/genes_cell2struct ${IMfilename}
fi

if [ "$ANNO_FORMAT" == '1' ]
then
    echo load the genome annotation in AGS format
	IMfilename=${ANNO_INPUT%.dat}.mat
    if [ ! -f ${IMfilename} ]
	then
		ln -s ${ANNO_INPUT%.dat}_files/genes.mat ${IMfilename}
		${DIR}/../bin/genes_cell2struct ${IMfilename}
    fi
fi

echo
echo %%%%%%%%%%%%%%%%%%%%
echo % 2. Read counting %
echo %%%%%%%%%%%%%%%%%%%%
echo

echo counting reads everlapping exons using given alignments

for REPLICATE_GROUP in $@
do
    IFS=':' read -a REPLICATES <<< "$REPLICATE_GROUP"
    for BAM_FILE in ${REPLICATES[@]}
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

tmpfile=`mktemp --tmpdir=/home/galaxy/tmp`

#${DIR}/../bin/get_read_counts ${GENES_FN} $tmpfile $@
echo "${DIR}/../bin/get_read_counts ${IMfilename} $tmpfile $@"
${DIR}/../bin/get_read_counts ${IMfilename} $tmpfile $@

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo testing genes for differential expression using given alignments

echo "cat ${DIR}/../src/difftest_deseq.R | $R_PATH --slave --args $tmpfile ${DESEQ_RES_FILE} $#"
cat ${DIR}/../src/difftest_deseq.R | $R_PATH --slave --args $tmpfile ${DESEQ_RES_FILE} $# 2> /dev/null

rm $tmpfile
rm ${tmpfile}_COUNTS.tab
rm ${tmpfile}_CONDITIONS.tab
echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
