#!/bin/sh
## Find all eucaryotes in the data dir and convert them
## using matlab script prepare_data.
##
## Works.

datadir="../../out/"
dirs=$(find $datadir -maxdepth 1 -mindepth 1 -type d \( ! -name .svn \))

for dir in $dirs
do
	# last directory is name of the organism
	organism=${dir##$datadir}
	echo "organism='$organism'"

	# call matlab prepare script
	matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir$organism'), exit;"
done

echo "Finished."
