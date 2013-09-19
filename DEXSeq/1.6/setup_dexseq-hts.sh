#!/bin/bash
set -e 

DIR=`dirname $0`
. ${DIR}/./bin/dexseq_config.sh

echo ==========================================
echo  DEXSeq-hts setup script \(DEXSeq version $DEXSEQ_VERSION\) 
echo ==========================================
echo
echo SAMTools directory \(currently set to \"$SAMTOOLS_DIR\", system version used if left empty\)
read SAMTOOLS_DIR
if [ "$SAMTOOLS_DIR" == "" ];
then
	if [ "$(which samtools)" != "" ] ;
	then
		SAMTOOLS_DIR=$(dirname $(which samtools)) 
	else
		echo samtools not found
		exit -1 ;
	fi
fi
echo '=>' Setting SAMTools directory to \"$SAMTOOLS_DIR\"
echo

echo Path to the python binary \(currently set to \"$PYTHON_PATH\", system version used, if left empty\)
read PYTHON_PATH
if [ "$PYTHON_PATH" == "" ];
then
    PYTHON_PATH=`which python`
	if [ "$PYTHON_PATH" == "" ];
	then
		echo python not found
		exit -1 
	fi
fi
echo '=>' Setting Python path to \"$PYTHON_PATH\"
echo

echo Path to HTSeq library files \(currently set to \"$PYTHONPATH\", system version is used if left empty\)
read PYTHONPATH
echo '=>' Setting HTSeq path to \"$PYTHONPATH\"
echo

echo Path to the R binary \(currently set to \"$R_PATH\", system version used, if left empty\)
read R_PATH
if [ "$R_PATH" == "" ];
then
    R_PATH=`which R`
	if [ "$R_PATH" == "" ];
	then
		echo R not found
		exit -1 
	fi
fi
echo '=>' Setting R path to \"$R_PATH\"
echo

cp -p bin/dexseq_config.sh bin/dexseq_config.sh.bk
grep -v -e SAMTOOLS_DIR -e PYTHON_PATH -e PYTHONPATH -e R_PATH -e $DEXSEQ_VERSION bin/dexseq_config.sh.bk > bin/dexseq_config.sh
echo
echo
echo generating config file

echo export DEXSEQ_VERSION=$DEXSEQ_VERSION >> bin/dexseq_config.sh
echo export SAMTOOLS_DIR=$SAMTOOLS_DIR >> bin/dexseq_config.sh
echo export PYTHON_PATH=$PYTHON_PATH >> bin/dexseq_config.sh
echo export PYTHONPATH=$PYTHONPATH >> bin/dexseq_config.sh
echo export R_PATH=$R_PATH >> bin/dexseq_config.sh

echo
echo Done.
echo 
