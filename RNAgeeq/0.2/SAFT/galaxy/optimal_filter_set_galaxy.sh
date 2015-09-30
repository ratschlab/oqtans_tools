#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    annotation="$1"
fi

if [ -z "$2" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    alignment="$2"
fi

if [ -z "$3" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    best_score="$3"
fi

if [ -z "$4" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    score_matrix="$4"
fi

if [ -z "$5" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in BAM> <best_score> <score_matrix> <filtered outfile>"
    echo ""
    exit 1
else
    filter_out="$5"
fi
SRC_DIR="$OQTANS_PATH/RNAgeeq/0.2/SAFT"
SAMTOOLS="$OQTANS_DEP_PATH/bin/"
export PYTHONPATH="$SRC_DIR/tools/:$PYTHONPATH"

echo "Generating annotation intron list"
echo ""

#if [ ! -f ${annotation}.introns ]
#then
    python $SRC_DIR/gen_intronlist_from_annotation.py -v -a $annotation -o ${annotation}.introns
#fi

echo "done"
echo ""


echo "Generating alignment feature list"
echo ""

#if [ ! -f ${alignment}.features ]
#then
    python $SRC_DIR/get_intron_features.py -b -a $alignment -o ${alignment}.features -s $SAMTOOLS
#fi

echo "done"
echo ""

echo "Search for optimal filter setting"
echo ""
python  $SRC_DIR/find_optimal_param_set.py -v -i ${annotation}.introns -f ${alignment}.features -b $best_score -m $score_matrix
echo "done"
echo ""

if [ ! -z "$filter_out" ]
then
    echo "Filter Alignment"
    echo ""
    min_ex_len=`tail -n 1 $best_score | cut -f 1`
    max_mm=`tail -n 1 $best_score | cut -f 2`
    min_support=`tail -n 1 $best_score | cut -f 3`
    support_string=""

    if [ "$min_support" != "1" ]
    then
        python $SRC_DIR/filter_features.py -e $min_ex_len -X $max_mm -m $min_support -i ${alignment}.features -o ${alignment}.features_filtered
        support_string="-i ${alignment}.features_filtered"
    fi

    (python $SRC_DIR/filter_alignment.py -b -a $alignment -o $filter_out -e $min_ex_len -X $max_mm $support_string -s $SAMTOOLS 2>&1 ||(echo Writing BAM file failed 1>&2))

    echo done
    echo ""
fi
