#/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2011 Max Planck Society
#

export DESEQ_VERSION=1.6.1
export DESEQ_PATH=$OQTANS_PATH/DESeq/1.6/
export DESEQ_SRC_PATH=$DESEQ_PATH/src
export DESEQ_BIN_PATH=$DESEQ_PATH/bin
export INTERPRETER=octave
export MATLAB_BIN_PATH=
export MATLAB_MEX_PATH=
export MATLAB_INCLUDE_DIR=
export OCTAVE_BIN_PATH=${OQTANS_DEP_PATH}/octave/bin/octave
export OCTAVE_MKOCT=${OQTANS_DEP_PATH}/octave/bin/mkoctfile
export SAMTOOLS_DIR=${OQTANS_DEP_PATH}/bin/
export PYTHON_PATH=/usr/bin/python${OQTANS_PYTHON_VERSION}
export SCIPY_PATH=${OQTANS_DEP_PATH}/lib/python${OQTANS_PYTHON_VERSION}/site-packages/
export LD_LIBRARY_PATH=/lib:$LD_LIBRARY_PATH
export R_PATH=${OQTANS_DEP_PATH}/R/bin/R
export ENVIRONMENT=galaxy
