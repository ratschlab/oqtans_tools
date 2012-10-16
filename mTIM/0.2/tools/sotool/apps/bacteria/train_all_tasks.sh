#!/bin/sh
## Find all eucaryotes in the data dir and convert them
## using matlab script prepare_data.
##
## Works.

datadir="../../data"
dirs=$(find $datadir -maxdepth 1 -mindepth 1 -type d \( ! -name .svn \))

out_dir="../out/All_Tasks"

num_train="10"
num_vald="50"

params_c="[0.1 1 10]"

organism=""

for dir in $dirs
do  
    dir2=${dir##$datadir}
    organism="$organism,'$dir2'"
done

organism=${organism#?}
organism="{ $organism }"

# call matlab prepare script
echo --------------------------------------
parameters="'$datadir',$organism,'$out_dir',$num_train,$num_vald,$params_c,[]"
echo $parameters

matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"


echo "Finished."
