#/bin/bash
##
# Galaxy wrapper script for rDiff version 0.3
# Copyright (C) 2013 cBio Department Memorial Sloan-Kettering Cancer Center
#
# This program is free software; you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation; either version 3 of the License, or 
# (at your option) any later version.
##

set -e 

PROG=`basename $0`
DIR=`dirname $0`
. ${DIR}/../bin/rdiff_config.sh

echo
echo ${PROG}: Oqtans http://galaxy.cbio.mskcc.org Galaxy wrapper for rDiff version $RDIFF_VERSION.
echo
echo rDiff performs differential expression testing from RNA-Seq measurements.
echo 

if [ -z "$1" -o "$1" == '--help' ];
then
  echo Usage: $0 poisson\|param\|nonparam in.gff readLength out.resultfile out_result_dir S1_in.bam?:S2_in.bam?:   
  echo "     :" $0 poisson\|param\|nonparam in.gff readLength out.resultfile out_result_dir S1_R1_in.bam?S1_R2_in.bam?:S2_R1_in.bam?S2_R2_in.bam?:   
  echo "   or:" $0 --help
  echo 
  false
fi

if [ "$1" != 'poisson' -a "$1" != 'param' -a "$1" != 'nonparam' ];
then
  echo invalid parameter: $1
  echo 
  echo "For usage:" $0 --help
  false
fi

TEST_METH=$1 ## Test method 
shift
GFF_INPUT=$1 ## Genome annotation in GFF format 
shift
readlength=$1 ## Sequencing read length 
shift 
RDIFF_OUT=$1 ## rDiff result file
shift
RDIFF_OUT_PATH=$1 ## temp session working directory 
shift

if [ -d $RDIFF_OUT_PATH ]
then       
	echo "Using the extra file path as : $RDIFF_OUT_PATH"
else
	mkdir -p ${RDIFF_OUT_PATH} ## create the temp working directory
fi 

## changing tool working directory to the session path
cd ${RDIFF_OUT_PATH}

## Seperating the files according to the sample and correspondinf replicates. 
SAMPLE_LIST=()
for SAMPLES in $@
do
    IFS='?' 
    for BAM_FILE in ${SAMPLES}
    do
        if [ $BAM_FILE = ":" ]; ## samples are seperated with ':'
        then
            SAMPLE_LIST+=(${SAMPLE_FNAME%?}) ## samples are separating 
            SAMPLE_FNAME=""
            continue
        fi
        if [ ! -f ${BAM_FILE}.bai ]
        then
            echo "Indexing $BAM_FILE"
            ${RDIFF_SAMTOOLS_INCLUDE_DIR}/samtools index $BAM_FILE
        else
            echo "$BAM_FILE already indexed"
        fi
        La_fn=$BAM_FILE
        SAMPLE_FNAME="$SAMPLE_FNAME$BAM_FILE," ## adding a ',' between each BAM files.
    done
done

## rDiff execution call 
${DIR}/../bin/rdiff -o ${RDIFF_OUT_PATH} -d / -g ${GFF_INPUT} -a ${SAMPLE_LIST[0]} -b ${SAMPLE_LIST[1]} -m ${TEST_METH} -L ${readlength}

## rdiff out file
ln -fs ${RDIFF_OUT_PATH}/P_values_rDiff_"${TEST_METH/param/parametric}".tab ${RDIFF_OUT}
