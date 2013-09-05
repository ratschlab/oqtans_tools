#!/bin/bash
set -e 

DIR=`dirname $0`
. ${DIR}/./bin/deseq_config.sh

echo ==========================================
echo  DESeq-hts setup script \(DESeq version $DESEQ_VERSION\) 
echo ==========================================
echo
echo DESeq-hts base directory \(currently set to \"$DESEQ_PATH\", suggest to set to \"`pwd`\", used if left empty\)
read DESEQ_PATH
if [ "$DESEQ_PATH" == "" ];
then
	DESEQ_PATH=`pwd`
fi
echo '=>' Setting DESeq-hts base directory to \"$DESEQ_PATH\"
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

echo Path to Scipy library files \(currently set to \"$SCIPY_PATH\", system version is used if left empty\)
read SCIPY_PATH
echo '=>' Setting Scipy path to \"$SCIPY_PATH\"
echo

echo Which interpreter should be used \(\"octave\" or \"matlab\"\)
read INTERPRETER  
if [ "$INTERPRETER" != 'octave' -a  "$INTERPRETER" != 'matlab' ];
then
	echo Unrecognized choice: \"$INTERPRETER\"
	echo Aborting
	false
fi
echo '=>' Setting interpreter to \"$INTERPRETER\"
echo

if [ "$INTERPRETER" == 'octave' ];
then
	echo Please enter the full path to octave \(currently set to \"$OCTAVE_BIN_PATH\", system version used, if left empty\)
	read OCTAVE_BIN_PATH
	if [ "$OCTAVE_BIN_PATH" == "" ];
	then
	    OCTAVE_BIN_PATH=`which octave` 
		if [ "$OCTAVE_BIN_PATH" == "" ];
		then
			echo octave not found
			exit -1
		fi
	fi
	echo '=>' Setting octave\'s path to \"$OCTAVE_BIN_PATH\"
	echo
	echo Please enter the full path to mkoctfile \(currently set to \"$OCTAVE_MKOCT\", system version used, if left empty\)
	read OCTAVE_MKOCT
	if [ "$OCTAVE_MKOCT" == "" ];
	then
	    OCTAVE_MKOCT=`which mkoctfile` 
		if [ "$OCTAVE_MKOCT" == "" ];
		then
			OCTAVE_MKOCT=$(dirname $OCTAVE_BIN_PATH)/mkoctfile
			if [ ! -f OCTAVE_MKOCT ];
			then
				echo mkoctfile not found
				exit -1
			fi
		fi
	fi
	echo '=>' Setting mkoctfile\'s path to \"$OCTAVE_MKOCT\"
	echo
	MATLAB_BIN_PATH=
fi
if [ "$INTERPRETER" == 'matlab' ];
then
	echo Please enter the full path to matlab \(currently set to \"$MATLAB_BIN_PATH\", system version used, if left empty\)
	read MATLAB_BIN_PATH
	if [ "${MATLAB_BIN_PATH}" == "" ];
	then
		MATLAB_BIN_PATH=`which matlab`
		if [ "$MATLAB_BIN_PATH" == "" ];
		then
			echo matlab not found
			exit -1
		fi
	fi
	if [ ! -f $MATLAB_BIN_PATH ];
	then
		echo matlab not found
		exit -1
	fi
	echo '=>' Setting matlab\'s path to \"$MATLAB_BIN_PATH\"
	echo
	echo Please enter the full path to mex binary \(currently set to \"$MATLAB_MEX_PATH\", system version used if left empty\)
	read MATLAB_MEX_PATH
	if [ "$MATLAB_MEX_PATH" == "" ];
	then
		MATLAB_MEX_PATH=`which mex`
		if [ "$MATLAB_MEX_PATH" == "" ];
		then
			echo mex not found
			exit -1
		fi
	fi
	if [ ! -f "$MATLAB_MEX_PATH" ];
	then
		echo mex not found
		exit -1
	fi
	echo '=>' Setting mex\' path to \"$MATLAB_MEX_PATH\"
	echo
	echo Please enter the full path to the matlab include directory \(currently set to \"$MATLAB_INCLUDE_DIR\", system version used, if left empty\)
	read MATLAB_INCLUDE_DIR
	if [ "$MATLAB_INCLUDE_DIR" == "" ];
	then
		MATLAB_INCLUDE_DIR=$(dirname $MATLAB_BIN_PATH)/../extern/include
	fi
	if [ ! -d "$MATLAB_INCLUDE_DIR" ];
	then
		echo matlab include dir not found
		exit -1
	fi
	echo '=>' Setting matlab\'s include directory to \"$MATLAB_INCLUDE_DIR\"
	echo
	OCTAVE_BIN_PATH=
fi

cp -p bin/deseq_config.sh bin/deseq_config.sh.bk
grep -v -e OCTAVE_BIN_PATH -e OCTAVE_MKOCT -e MATLAB_BIN_PATH -e MATLAB_MEX_PATH -e MATLAB_INCLUDE_DIR \
    -e DESEQ_PATH -e DESEQ_SRC_PATH -e DESEQ_BIN_PATH \
    -e INTERPRETER -e SAMTOOLS_DIR -e PYTHON_PATH -e SCIPY_PATH -e R_PATH -e $DESEQ_VERSION bin/deseq_config.sh.bk  \
    > bin/deseq_config.sh
echo
echo
echo generating config file

echo export DESEQ_VERSION=$DESEQ_VERSION >> bin/deseq_config.sh
echo export DESEQ_PATH=$DESEQ_PATH >> bin/deseq_config.sh
echo export DESEQ_SRC_PATH=${DESEQ_PATH}/src >> bin/deseq_config.sh
echo export DESEQ_BIN_PATH=${DESEQ_PATH}/bin >> bin/deseq_config.sh
echo export INTERPRETER=$INTERPRETER >> bin/deseq_config.sh
echo export MATLAB_BIN_PATH=$MATLAB_BIN_PATH >> bin/deseq_config.sh
echo export MATLAB_MEX_PATH=$MATLAB_MEX_PATH >> bin/deseq_config.sh
echo export MATLAB_INCLUDE_DIR=$MATLAB_INCLUDE_DIR >> bin/deseq_config.sh
echo export OCTAVE_BIN_PATH=$OCTAVE_BIN_PATH >> bin/deseq_config.sh
echo export OCTAVE_MKOCT=$OCTAVE_MKOCT >> bin/deseq_config.sh
echo export SAMTOOLS_DIR=$SAMTOOLS_DIR >> bin/deseq_config.sh
echo export PYTHON_PATH=$PYTHON_PATH >> bin/deseq_config.sh
echo export SCIPY_PATH=$SCIPY_PATH >> bin/deseq_config.sh
echo export R_PATH=$R_PATH >> bin/deseq_config.sh

echo
echo Done.
echo 
