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



export DESEQ_VERSION=1.6.0
export DESEQ_PATH=/home/galaxy/galaxy-2.1.2009/tools/quantitation_tools/deseq
export DESEQ_SRC_PATH=$DESEQ_PATH/src
export DESEQ_BIN_PATH=$DESEQ_PATH/bin
export INTERPRETER=octave
export MATLAB_BIN_PATH=/home/galaxy/software/matlab-7.6/bin/matlab
export MATLAB_MEX_PATH=/home/galaxy/software/matlab-7.6/bin/mex
export MATLAB_INCLUDE_DIR=/home/galaxy/software/matlab-7.6/extern/include
export OCTAVE_BIN_PATH=/home/galaxy/software/octave-3.2.3-64-new/bin/octave
export OCTAVE_MKOCT=/home/galaxy/software/octave-3.2.3-64-new/bin/mkoctfile
export SAMTOOLS_DIR=/home/galaxy/software/samtools
#export PYTHON_PATH=/home/galaxy/galaxy-2.1.2009/galaxy-python/python
export PYTHON_PATH=/usr/bin/python
export SCIPY_PATH=/home/galaxy/python/lib/python2.5/site-packages/
export LD_LIBRARY_PATH=/home/galaxy/software/lib:/home/galaxy/lib:$LD_LIBRARY_PATH
#export R_PATH=/usr/bin/R
export R_PATH=/home/galaxy/software/R-2.14.1/bin/R
export ENVIRONMENT=galaxy
