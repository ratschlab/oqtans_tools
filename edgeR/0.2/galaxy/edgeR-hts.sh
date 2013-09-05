#/bin/bash
#
# Under development 
#

set -e 

PROG=`basename $0`
DIR=`dirname $0`

echo
echo
echo edgeR performs differential expression testing from RNA-Seq measurements.
echo 

ANNO_INPUT=${1}
shift
EDGER_RES_FILE=${1}
shift
GENES_FN=${1}
shift

mkdir -p `dirname $GENES_FN`

echo %%%%%%%%%%%%%%%%%%%%%%%
echo % 1. Data preparation %
echo %%%%%%%%%%%%%%%%%%%%%%%
echo
echo load the genome annotation in GFF3 format and create an annotation object
echo
export PYTHONPATH=$PYTHONPATH:${SCIPY_PATH}
echo ${ANNO_INPUT}
${PYTHON_PATH} ${DIR}/../tools/ParseGFF.py ${ANNO_INPUT} ${GENES_FN}
${DIR}/../bin/genes_cell2struct ${GENES_FN} 2>&1
echo 
echo genome annotation stored in $GENES_FN

echo
echo %%%%%%%%%%%%%%%%%%%%
echo % 2. Read counting %
echo %%%%%%%%%%%%%%%%%%%%
echo

echo counting reads overlapping exons using given alignments
for REPLICATE_GROUP in $@
do
    IFS=':'
    for BAM_FILE in ${REPLICATE_GROUP}
    do
        echo
        if [ ! -f ${BAM_FILE}.bai ]
        then
            echo "Indexing $BAM_FILE"
            ${SAMTOOLS_DIR}/samtools index $BAM_FILE
        else
            echo "$BAM_FILE already indexed"
        fi
        echo
    done
done
tmpfile=`mktemp --tmpdir=/tmp`

echo "${DIR}/../bin/get_read_counts ${GENES_FN} $tmpfile $@"
${DIR}/../bin/get_read_counts ${GENES_FN} $tmpfile $@ 2>&1

echo $tmpfile ${tmpfile}_COUNTS.tab ${tmpfile}_CONDITIONS.tab

echo
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo % 3. Differential testing %
echo %%%%%%%%%%%%%%%%%%%%%%%%%%%
echo
#
#echo testing genes for differential expression using given alignments
#
echo "cat ${DIR}/../src/edger-hts.R | $R_PATH --slave --args $tmpfile ${DESEQ_RES_FILE} $#"
cat ${DIR}/../src/edger-hts.R| $R_PATH --slave --args $tmpfile ${DESEQ_RES_FILE} $# 2> /dev/null
#
##rm $tmpfile ${tmpfile}_COUNTS.tab ${tmpfile}_CONDITIONS.tab
#echo $tmpfile ${tmpfile}_COUNTS.tab ${tmpfile}_CONDITIONS.tab
#
echo
echo %%%%%%%%
echo % Done %
echo %%%%%%%%
echo
