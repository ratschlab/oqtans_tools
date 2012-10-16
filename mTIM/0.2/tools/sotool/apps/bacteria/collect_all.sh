#!/bin/sh
##
## Works.

datadir="../../out/$1"
dirs=$(find $datadir -mindepth 1 -iname 'evaluation.mat' -type f \( ! -name .svn \))

organism=""

for dir in $dirs
do  
    organism="$organism,'$dir'"
done

organism=${organism#?}
organism="{ $organism }"

# call matlab prepare script
echo --------------------------------------
echo $organism
matlab -nojvm -nosplash -nodesktop -r "collect_evals($organism), exit;"


echo "Finished."
