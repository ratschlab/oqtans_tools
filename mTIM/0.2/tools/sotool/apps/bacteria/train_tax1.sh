#!/bin/sh
##
## Trains and evaluates the mini 4-taxonomy.
## Assumes root node is already trained and tested.
##

num_train="10"
num_vald="20"


datadir="../../data"

out_dir_root="../out/tax1_$num_train/root"

out_dir_b="../out/tax1_$num_train/bacillus"
out_dir_e="../out/tax1_$num_train/escherichia"

out_dir_ba="../out/tax1_$num_train/b.anthracis"
out_dir_bs="../out/tax1_$num_train/b.subtilis"
out_dir_ec="../out/tax1_$num_train/e.coli"
out_dir_ef="../out/tax1_$num_train/e.fergusonii"



params_c="[0.1 1 10]"
params_b="[0.01 0.1 0.5 0.8 1.0]"

name_ba="'Bacillus_anthracis_Ames_uid57909'"
name_bs="'Bacillus_subtilis_168_uid57675'"
name_ec="'Escherichia_coli_BW2952_uid59391'"
name_ef="'Escherichia_fergusonii_ATCC_35469_uid59375'"

node0="{$name_ba, $name_bs, $name_ec, $name_ef}"

node11="{$name_ba, $name_bs}"
node12="{$name_ec, $name_ef}"

node21="{$name_ba}"
node22="{$name_bs}"
node23="{$name_ec}"
node24="{$name_ef}"

root_file_11="$out_dir_b/evaluation.mat"
root_file_12="$out_dir_e/evaluation.mat"

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
parameters="'$datadir',$node11,'$out_dir_b',$num_train,$num_vald,$params_c,$params_b,'../$root_file'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_b'), exit;"


parameters="'$datadir',$node12,'$out_dir_e',$num_train,$num_vald,$params_c,$params_b,'../$root_file'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_e'), exit;"




echo --------------------------------------
parameters="'$datadir',$node21,'$out_dir_ba',$num_train,$num_vald,$params_c,$params_b,'../$root_file_11'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_ba'), exit;"


parameters="'$datadir',$node22,'$out_dir_bs',$num_train,$num_vald,$params_c,$params_b,'../$root_file_11'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_bs'), exit;"


parameters="'$datadir',$node23,'$out_dir_ec',$num_train,$num_vald,$params_c,$params_b,'../$root_file_12'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_ec'), exit;"

parameters="'$datadir',$node24,'$out_dir_ef',$num_train,$num_vald,$params_c,$params_b,'../$root_file_12'"
matlab -nojvm -nosplash -nodesktop -r "train_orgs($parameters), exit;"
echo "EVALUATION"
matlab -nojvm -nosplash -nodesktop -r "evaluate_result('$datadir/$out_dir_ef'), exit;"


echo "Finished."
