#!/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
#

set -e 

. ./bin/rquant_config.sh

echo =====================================
echo  rQuant setup script \(version $RQUANT_VERSION\) 
echo =====================================
echo

echo rQuant base directory \(currently set to \"$RQUANT_PATH\", suggest to set to \"`pwd`\", used if left empty\)
read RQUANT_PATH       
if [ "$RQUANT_PATH" == "" ];
then
	RQUANT_PATH=`pwd`
fi
echo '=>' Setting rQuant base directory to \"$RQUANT_PATH\"
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

echo Path to Scipy installation \(currently set to \"$SCIPY_PATH\", system version is used if left empty\)
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
	echo '=>' Setting octave\'s path to \"$OCTAVE_MKOCT\"
	echo

	MATLAB_BIN_PATH=
fi

if [ "$INTERPRETER" == 'matlab' ];
then
	echo Please enter the full path to matlab \(currently set to \"$MATLAB_BIN_PATH\", system version used, if left empty\)
	read MATLAB_BIN_PATH
	if [ "$MATLAB_BIN_PATH" == "" ];
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

cp -p bin/rquant_config.sh bin/rquant_config.sh.bak
grep -v -e OCTAVE_BIN_PATH -e OCTAVE_MKOCT -e MATLAB_BIN_PATH -e MATLAB_MEX_PATH -e MATLAB_INCLUDE_DIR \
    -e RQUANT_PATH -e RQUANT_SRC_PATH -e RQUANT_BIN_PATH \
    -e INTERPRETER -e SAMTOOLS_DIR -e PYTHON_PATH -e SCIPY_PATH -e $RQUANT_VERSION bin/rquant_config.sh.bak  \
    > bin/rquant_config.sh

echo
echo
echo generating config file

# appending the relevant lines to rquant_config.sh
echo export RQUANT_VERSION=$RQUANT_VERSION >> bin/rquant_config.sh
echo export RQUANT_PATH=$RQUANT_PATH >> bin/rquant_config.sh
echo export RQUANT_SRC_PATH=${RQUANT_PATH}/src >> bin/rquant_config.sh
echo export RQUANT_BIN_PATH=${RQUANT_PATH}/bin >> bin/rquant_config.sh
echo export INTERPRETER=$INTERPRETER >> bin/rquant_config.sh
echo export MATLAB_BIN_PATH=$MATLAB_BIN_PATH >> bin/rquant_config.sh
echo export MATLAB_MEX_PATH=$MATLAB_MEX_PATH >> bin/rquant_config.sh
echo export MATLAB_INCLUDE_DIR=$MATLAB_INCLUDE_DIR >> bin/rquant_config.sh
echo export OCTAVE_BIN_PATH=$OCTAVE_BIN_PATH >> bin/rquant_config.sh
echo export OCTAVE_MKOCT=$OCTAVE_MKOCT >> bin/rquant_config.sh
echo export SAMTOOLS_DIR=$SAMTOOLS_DIR >> bin/rquant_config.sh
echo export PYTHON_PATH=$PYTHON_PATH >> bin/rquant_config.sh
echo export SCIPY_PATH=$SCIPY_PATH >> bin/rquant_config.sh

echo
echo
echo compiling source code
cd mex
if [ "$INTERPRETER" == "octave" ];
then
	make octave
else
	make matlab
fi
cd ..

echo
echo Done.
echo 
