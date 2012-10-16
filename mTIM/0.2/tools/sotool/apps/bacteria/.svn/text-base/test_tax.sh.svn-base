#!/bin/sh
##
## Trains and evaluates the mini 4-taxonomy.
## Assumes root node is already trained and tested.
##

num_train="10"
num_vald="20"


datadir="../../data"

out_dir_root="../out/tax1_$num_train/root"

out_dir_bs="../out/tax1_$num_train/b.subtilis"
out_dir_ec="../out/tax1_$num_train/e.coli"



params_c="[1 10]"
params_b="[0.01 1.0]"

name_bs="'Bacillus_subtilis_168_uid57675'"
name_ec="'Escherichia_coli_BW2952_uid59391'"

node0="{$name_bs, $name_ec}"

node11="{$name_bs}"
node12="{$name_ec}"

root_file="$out_dir_root/evaluation.mat"

## call matlab scripts for training and evaluation
##
##
echo --------------------------------------
parameters="'$datadir',$node0,'$out_dir_root',$num_train,$num_vald,$params_c,[]"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_root'), exit;"



echo --------------------------------------
parameters="'$datadir',$node11,'$out_dir_bs',$num_train,$num_vald,$params_c,$params_b,'../$root_file'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_bs'), exit;"


parameters="'$datadir',$node12,'$out_dir_ec',$num_train,$num_vald,$params_c,$params_b,'../$root_file'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_ec'), exit;"


echo "Finished."
