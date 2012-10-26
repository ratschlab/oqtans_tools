#!/bin/bash
#	Main wrapper script to run mTIM program 

# 	mTIM_Train requires 
# 	read alignment file in BAM 
# 	genome sequence in FASTA
#	acc/don splice site prediction in SPF format 
# 	partial genome annotation in GFF3 
set -e 
DIR=`dirname $0`

FASTA_IN=${1}
GFF_IN=${2}
BAM_IN=${3}
ACC_SPF=${4}
DON_SPF=${5}
#OUT_RES_FILE=${6}
OUT_DIR=${6}
PROG='mTIM-Train'
echo
echo ${PROG}: This program is part of the mTIM version 0.2
echo
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '1. Genome Sequence Preparation'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
/usr/bin/python ${DIR}/../tools/GenomePrep.py $FASTA_IN $OUT_DIR 
echo '  Done...'
echo
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '2. Genome Annotation Preparation'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
/usr/bin/python ${DIR}/../tools/ParseGFF.py $GFF_IN $OUT_DIR/genes.mat
${DIR}/../bin/genes_cell2struct ${OUT_DIR}/
echo '  Done...'
echo
echo '%%%%%%%%%%%'
echo '3. Training'
echo '%%%%%%%%%%%'
### check if bam file is indexed
if [ ! -f ${BAM_IN}.bai ]
then
	$OQTANS_DEP_PATH/bin/samtools index ${BAM_IN}
fi
${DIR}/../bin/mTIM_galaxy_train ${OUT_DIR}/ ${ACC_SPF} ${DON_SPF} ${BAM_IN} ${OUT_DIR}/ ${OUT_DIR}/ 
echo '  Done...'
echo
