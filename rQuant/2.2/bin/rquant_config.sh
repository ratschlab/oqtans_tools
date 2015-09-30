#!/bin/bash
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
export RQUANT_VERSION=2.2
export RQUANT_PATH=${OQTANS_PATH}/rQuant/2.2
export RQUANT_SRC_PATH=${RQUANT_PATH}/src
export RQUANT_BIN_PATH=${RQUANT_PATH}/bin
export INTERPRETER=octave
export MATLAB_BIN_PATH=
export MATLAB_MEX_PATH=
export MATLAB_INCLUDE_DIR=
export OCTAVE_BIN_PATH=${OQTANS_DEP_PATH}/octave/bin/octave
export OCTAVE_MKOCT=${OQTANS_DEP_PATH}/octave/bin/mkoctfile
export SAMTOOLS_DIR=${OQTANS_DEP_PATH}/bin/
export PYTHON_PATH=${OQTANS_PYTHON}
export SCIPY_PATH=${OQTANS_DEP_PATH}/lib/python${OQTANS_PYTHON_VERSION}/site-packages/:${OQTANS_DEP_PATH}/lib64/python${OQTANS_PYTHON_VERSION}/site-packages/:$PYTHONPATH
