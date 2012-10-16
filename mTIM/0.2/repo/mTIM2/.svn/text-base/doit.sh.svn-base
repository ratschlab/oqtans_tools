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

echo "Experiment name (and prefix of config file) is ${prefix}"
echo "Copy src directory to ${dest}"

## create directory
mkdir ~/${dest}
cp -r -L src ~/${dest}

##echo "Submitting job."
##qsub -p 999 remote.sh ${prefix} ${exp}

