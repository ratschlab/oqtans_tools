#!/bin/bash
#
# oqtans submodule compile script
# 
set -e

basedir=`dirname $0`
if [ "$basedir" == "." ];
then
    basedir=`pwd`/
fi

echo $basedir

for i in PALMapper/0.5 DESeq/1.12 DESeq2/1.0.19 rQuant/2.2 mTIM/0.2 Trinity/r2013_08_14 BWA/0.6.2 edgeR/0.2 EasySVM/0.3.3 rDiff/0.3
do
    cd $basedir/$i
    echo ==============================================================
    echo ===================== making $i ===============
    echo ==============================================================
    ./oqtans_make.sh $1
    echo ================== created binary files for $i ===============
done
