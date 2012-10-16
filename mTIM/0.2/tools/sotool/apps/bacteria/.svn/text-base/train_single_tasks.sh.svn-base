#!/bin/sh
## Find all eucaryotes in the data dir and convert them
## using matlab script prepare_data.
##
## Works.

datadir="../../data/"
dirs=$(find $datadir -maxdepth 1 -mindepth 1 -type d \( ! -name .svn \))

echo $dirs

num_train="40";
num_vald="120";
params_c="[1 10 100 250]"

for dir in $dirs
do
	# last directory is name of the organism
	organism=${dir##$datadir}
	##echo "organism='$organism'"

	## check if directories already exist
	if [ ! -d $datadir/../out/${organism}_$num_train ]
	then
	    ## doesnt exist -> train
	    echo "$datadir/../out/${organism}_$num_train"
	    # call matlab prepare script
	    matlab -nojvm -nosplash -nodesktop -r "train_orgs('$datadir',{'$organism'},'../out/${organism}_$num_train',$num_train,$num_vald,$params_c,[]), exit;"
	    echo "$datadir/../out/${organism}_$num_train"
	    matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/../out/${organism}_$num_train'), exit;"
        fi

	# call matlab prepare script
	#matlab -nojvm -nosplash -nodesktop -r "train_orgs('$datadir',{'$organism'},'../out/${organism}_$num_train',$num_train,$num_vald,$params_c,[]), exit;"
	#echo "$datadir/../out/${organism}_$num_train"
	#matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/../out/${organism}_$num_train'), exit;"
done

echo "Finished."
