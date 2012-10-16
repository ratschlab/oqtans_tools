#!/bin/bash

set -e

basedir=`dirname $0`
if [ "$basedir" == "." ];
then 
    basedir=`pwd`/
fi

echo $basedir

for i in palmapper/0.5 DESeq/1.6 rDiff/0.1 rQuant/2.2 mTIM/0.2 trinity/r2012-06-08 bwa/0.6.2
do
    cd $basedir/$i
    echo ==============================================================
    echo ===================== making $i ===============
    echo ==============================================================
    ./oqtans_make.sh $1
done


