#!/bin/bash

set -e

if [ -z "$1" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in SAM>"
    exit 1
else
    annotation="$1"
fi

if [ -z "$2" ]
then
    echo "usage: $0 <annotation in GFF3> <alignment in SAM>"
    exit 1
else
    alignment="$2"
fi


echo "Generating annotation intron list"
echo ""

if [ ! -f ${annotation}.introns ]
then
    python gen_intronlist_from_annotation.py -v -a $annotation -o ${annotation}.introns
fi

echo "done"
echo ""


echo "Gernerating alignment feature list"
echo ""

if [ ! -f ${alignment}.features ]
then
    python get_intron_features.py -v -a $alignment -o ${alignment}.features
fi

echo "done"
echo ""

echo "Evaluating alignment"
echo ""
python evaluate_features.py -v -i ${annotation}.introns -f ${alignment}.features
echo "done"
echo ""
