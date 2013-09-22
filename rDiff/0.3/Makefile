include bin/rdiff_config.sh
# Copyright (C) 2009-2013 Max Planck Society, Sloan-Kettering Institute

all:	mexfiles

mexfiles:
	echo Entering ./mex
	cd mex ; make octave
	echo Entering ./src/locfit/Source
	cd src/locfit/Source ; make octave

clean:	
	echo Entering ./mex
	cd mex ; make clean
	echo Entering ./examples
	cd examples ; make clean
	echo Entering ./src/locfit/Source
	cd src/locfit/Source ; make clean
	cd src/locfit/mex/ ; rm *

example:
	echo Entering ./examples
	cd examples ; make example

more_examples:	
	echo Entering ./examples
	cd examples ; make threeexamples
