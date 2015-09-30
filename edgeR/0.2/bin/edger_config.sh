#!/bin/bash
export edgeR_VERSION=3.2.4
export edgeR_PATH=${OQTANS_PATH}/edgeR/0.2
export edgeR_SRC_PATH=${edgeR_PATH}/src
export edgeR_BIN_PATH=${edgeR_PATH}/bin
export INTERPRETER=octave
export MATLAB_BIN_PATH=
export MATLAB_MEX_PATH=
export MATLAB_INCLUDE_DIR=
export OCTAVE_BIN_PATH=${OQTANS_DEP_PATH}/octave/bin/octave
export OCTAVE_MKOCT=${OQTANS_DEP_PATH}/octave/bin/mkoctfile
export SAMTOOLS_DIR=${OQTANS_SRC_PATH}/samtools-0.1.19
export PYTHON_PATH=${OQTANS_PYTHON}
export SCIPY_PATH=${OQTANS_DEP_PATH}/lib/python${OQTANS_PYTHON_VERSION}/site-package:${OQTANS_DEP_PATH}/lib64/python${OQTANS_PYTHON_VERSION}/site-package:$PYTHONPATH
export R_PATH=${OQTANS_R}
export LD_LIBRARY_PATH=$OQTANS_DEP_PATH/lib:$OQTANS_DEP_PATH/lib64::$LD_LIBRARY_PATH
export ENVIRONMENT=galaxy

