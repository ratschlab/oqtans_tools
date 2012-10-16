#!/bin/bash
# $1 samtoolsdr
# $2 tracks_dir
# $3 genome_fasta
# [$4] sam file
#
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2010 Regina Bohnert
# Copyright (C) 2009-2010 Max Planck Society
#


set -e

samtoolsdr=$1/

fastaname=$3

if [ -e ${fastaname}.fai ]
then
    echo "${fastaname}.fai already exists"
else
    echo "${fastaname}.fai is generated"
    $samtoolsdr./samtools faidx $fastaname
fi

if [ $# -lt 4 ]; then 
    sam_files=`ls ${2}/*.sam`
else
    sam_files=$4
fi
    
for n in $sam_files
do
    m=`echo $n | sed "s/.sam/.bam/g"`
    if [ ! -f $m ]
    then
		echo "creating file $m"
		$samtoolsdr./samtools import ${fastaname}.fai $n $m
		# BAM must be sorted by start position to use random access
		$samtoolsdr./samtools sort $m ${m}_sorted
		mv ${m}_sorted.bam $m
		# index for BAM file
		$samtoolsdr./samtools index $m
    fi
done
