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

if [ -z "$5" -o "$1" == '--help' ];
then
  echo Usage: $0
  echo "   or:" $0 --help
  false
fi 

ANNO_INPUT=${1}
shift
ANNO_FORMAT=${1}
shift
DESEQ_RES_FILE=${1}
shift


echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo

#Matlab structure file for genome annotation
ANNO_MAT="$(mktemp)"
mv $ANNO_MAT $ANNO_MAT.mat
ANNO_MAT=$ANNO_MAT.mat
echo $ANNO_MAT

#Annotation in GFF3 format
if [ "$ANNO_FORMAT" != "1" ]
then
    echo load the genome annotation in GFF3 format and create an annotation object
    echo
    export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
    ${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${ANNO_INPUT} all all $ANNO_MAT
    ${DIR}/../bin/genes_cell2struct $ANNO_MAT
fi

#Annotation in AGS format
if [ "$ANNO_FORMAT" == '1' ]
then
    echo load the genome annotation in AGS format
    ln -s ${ANNO_INPUT} ${ANNO_MAT}
    ${DIR}/../bin/genes_cell2struct ${ANNO_MAT}
fi

echo
echo %%%%%%%%%%%%%%%%%%%%
echo % 2. Read counting %
echo %%%%%%%%%%%%%%%%%%%%
echo

echo counting reads everlapping exons using given alignments

for BAM_FILE in $@
do	
	echo
	echo "Indexing bam file $BAM_FILE"
	echo
	if [ ! -f ${BAM_FILE}.bai ]
	then
		${SAMTOOLS_DIR}/samtools index $BAM_FILE
	fi
done

#tmpfile=`mktemp --tmpdir=/home/galaxy/tmp`
tmpfile=`tempfile -d /tmp`

echo "${DIR}/../bin/get_read_counts ${ANNO_MAT} $tmpfile $@"
${DIR}/../bin/get_read_counts ${ANNO_MAT} $tmpfile $@

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo

echo testing genes for differential expression using given alignments
echo "cat ${DIR}/../src/difftest_deseq.R \| /usr/bin/R --slave --args $tmpfile ${DESEQ_RES_FILE} $#"
cat ${DIR}/../src/difftest_deseq.R | /usr/bin/R --slave --args $tmpfile ${DESEQ_RES_FILE} $# 2> /dev/null

rm $tmpfile
rm ${tmpfile}_COUNTS.tab
echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
