#!/bin/bash
set -e
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${OQTANS_DEP_PATH}/octave/lib/
OCTAVE_BIN_PATH=${OQTANS_DEP_PATH}/octave/bin/octave
SRC_PATH=${OQTANS_PATH}/mTIM/0.2/src

echo exit | ${OCTAVE_BIN_PATH} --eval "global SHELL_INTERPRETER_INVOKE; SHELL_INTERPRETER_INVOKE=1; addpath $SRC_PATH; mTIM_config; $1($2); exit;" || (echo starting Octave failed; exit -1) ;
