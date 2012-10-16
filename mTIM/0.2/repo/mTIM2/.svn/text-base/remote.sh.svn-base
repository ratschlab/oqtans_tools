#!/bin/sh

## Copies the current source of mTIM to a specific
## loction defined in 'dest' and starts the whole
## pipeline of mtim using qsub command.
##
## Arg1 : config file prefix
## Arg2 : experiment name (destination output directory)
##
##
## written by nico goernitz 2012

prefix=$1;
exp=$2;
dest="mtim/${exp}";


/fml/ag-raetsch/share/software/matlabR2011a/bin/matlab -nodesktop -nosplash -nojvm -r "cd ~/${dest}/src; mTIM_all('$prefix'); exit;"


