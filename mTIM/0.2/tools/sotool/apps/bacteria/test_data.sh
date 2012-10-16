#!/bin/sh
## Find all eucaryotes in the data dir and do some simple
## checking whether the data is okay or there might be some
## errors.
##

datadir="../../data/"
dirs=$(find $datadir -maxdepth 1 -mindepth 1 -type d \( ! -name .svn \))

organism=""

for dir in $dirs
do  
    organism="$organism,'$dir'"
done

organism=${organism#?}
organism="{ $organism }"

# call matlab check_data script
echo --------------------------------------
echo $organism
matlab -nojvm -nosplash -nodesktop -r "check_data($organism), exit;"


echo "Finished."
