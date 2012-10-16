#!/bin/sh
## Find all eucaryotes in the data dir and convert them
## using matlab script prepare_data.
##
## Works.


datadir="../../data/"
dirs=$(find $datadir -maxdepth 1 -mindepth 1 -type d \( ! -name .svn \))


## Give some help
echo ---------------------------------------------------------------------------
echo Make sure that the data of the organisms lie in the directory: $datadir
echo And that there is a valid idx_to_codon.mat in the current directory. 
echo
echo The data has to be pre-processed and should have the following structure:
echo  $datadir/organism/mat/example_*.mat
echo Then results will be stored in:
echo  $datadir/organism/data.mat
echo
echo Possible Problems:
echo 1. There is no /mat -directory
echo Solution: Let the data be prepared by Christian Widmer
echo ---------------------------------------------------------------------------
echo 



for dir in $dirs
do
	# last directory is name of the organism
	organism=${dir##$datadir}
	echo "organism='$organism'"

	# call matlab prepare script
	matlab -nosplash -nojvm -nodesktop -r "prepare_data('$datadir','$organism'), exit;"
done

echo "Prepared: $dirs"
echo "Finished."
