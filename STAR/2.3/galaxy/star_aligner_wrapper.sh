#!/bin/bash
##
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Copyright (C) 2013 Memorial Sloan-Kettering Cancer Center 
##

## Can adjust this as appropriate for the system.
OPTIONS=" --genomeLoad NoSharedMemory"

## getting arguments from wrapper script 
GenomeSource=${1}
shift

if [ ${GenomeSource} = "history" ];
then 
    echo
    echo "STAR runMode generate Genome index"
    echo
    genomeFILE=${1}
    shift

    tmpDIR=`mktemp -d --tmpdir=/tmp`
    STAR --runMode genomeGenerate --genomeDir ${tmpDIR} --genomeFastaFiles ${genomeFILE} --runThreadN 4 
    wait
    OPTIONS=$OPTIONS" --genomeDir $tmpDIR"    
else
    ## built-in index  
    genomeDIR=${1}
    shift

    OPTIONS=$OPTIONS" --genomeDir $genomeDIR"    
fi 

LibTyp=${1}
shift 

FQ1=${1}
shift 

OPTIONS=$OPTIONS" --readFilesIn $FQ1"

if [ ${LibTyp} = "paired" ];
then
    FQ2=${1}
    shift
    OPTIONS=$OPTIONS" $FQ2"
fi 

AnnoFtype=${1}
shift 

ANNOFile=${1}
shift 

if [ ${AnnoFtype} = "bed" ];
then
    OPTIONS=$OPTIONS" --sjdbFileChrStartEnd $ANNOFile"
elif [ ${AnnoFtype} = "gtf" ]
then 
    OPTIONS=$OPTIONS" --sjdbGTFfile $ANNOFile"
else
    OPTIONS=$OPTIONS" --sjdbGTFfile $ANNOFile --sjdbGTFtagExonParentTranscript Parent"
fi

## parallel run settings 
OPTIONS=$OPTIONS" --runThreadN 2"

## BAM conversion at working path  
OUTPATH=${1}
shift 

mkdir -p $OUTPATH
cd $OUTPATH

mkfifo Aligned.out.sam
(cat Aligned.out.sam | samtools view -Shb - | samtools sort - Aligned.sort.out) 2>&1 &

defaultSET=${1}
shift 
if [ $defaultSET != "full" ];
then 
    echo "STAR run with default settings"
    STAR ${OPTIONS} || echo ERROR: STAR failed
else
    ## adding additional parameters 
    NM=${1}
    shift 
    OPTIONS=$OPTIONS" --outFilterMismatchNmax $NM"

    MMscore=${1}
    shift 
    OPTIONS=$OPTIONS" --outFilterMultimapScoreRange $MMscore"

    MM=${1}
    shift
    OPTIONS=$OPTIONS" --outFilterMultimapNmax $MM"

    maxIntron=${1}
    shift 
    OPTIONS=$OPTIONS" --alignIntronMax $maxIntron"

    mateGap=${1}
    shift 
    OPTIONS=$OPTIONS" --alignMatesGapMax $mateGap"

    sjScore=${1}
    shift
    OPTIONS=$OPTIONS" --sjdbScore $sjScore"

    sjOVH=${1}
    shift
    OPTIONS=$OPTIONS" --sjdbOverhang $sjOVH"

    #echo ${OPTIONS} 
    ## run STAR with the provided advanced settings 
    STAR ${OPTIONS} || echo ERROR: STAR failed
fi 
wait 

## moving file from working dir and cleaning up 
OUTFILE=${1}
shift
mv $OUTPATH/Aligned.sort.out.bam $OUTFILE

SPJNFL=${1}
shift
mv $OUTPATH/SJ.out.tab $SPJNFL

SUMMARY=${1}
shift
mv $OUTPATH/Log.final.out $SUMMARY

rm -fr $OUTPATH
