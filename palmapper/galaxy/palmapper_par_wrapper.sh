#!/bin/bash

set -e

date

export TEMP=/tmp
qsub=/opt/sge/bin/lx24-amd64/qsub
qstat=/opt/sge/bin/lx24-amd64/qstat

nlines=2000000

if [ "$1" == "-1" ]; 
then
    shift
    echo one input
    input1=${1#--input1=}
    if [ ! -f $input1 ];
    then
	echo error: file $input1 does not exist
	exit -1
    fi
    shift
else
    shift
    echo two inputs
    input1=${1#--input1=}
    if [ ! -f $input1 ];
    then
	echo error: file $input1 does not exist
	exit -1
    fi
    shift
    input2=${1#--input2=}
    if [ ! -f $input2 ];
    then
	echo error: file $input2 does not exist
	exit -2
    fi
    shift
fi

echo split -d -a 4 -l $nlines $input1 ${input1}.
split -d -a 4 -l $nlines $input1 ${input1}. 

if [ "$input2" != "" ];
then
    echo split -d -a 4 -l $nlines $input2 ${input2}.
    split -d -a 4 -l $nlines $input2 ${input2}.
fi

wait

suffix=`tempfile`
suffix=${suffix:(-5)}

basedir=`dirname $0`
for i in ${input1}.????; do
    splitno=${i:(-4)} 
    if [ "$input2" == "" ];
    then
	echo calling "palmapper on file $i"
	echo "./palmapper_wrapper.py --input1=$i --outsuffix=.${splitno}"
	echo "python $basedir/palmapper_wrapper.py --input1=$i --outsuffix=.${splitno} $* > /mnt/galaxyData/tmp/palmapper.$suffix.$splitno.log" | ${qsub} -cwd -N "plm$suffix$splitno" -pe "*" 1-4 -j y -o /mnt/galaxyData/tmp/palmapper.$suffix.$splitno.log
    else
	echo calling "./palmapper_wrapper.py --input1=${input1}.$splitno --input2=${input2}.$splitno "
	echo "python $basedir/palmapper_wrapper.py --input1=$i --input2=${input2}.$splitno --outsuffix=.${splitno} $* > /mnt/galaxyData/tmp/palmapper.$suffix.$splitno.log"
	echo "python $basedir/palmapper_wrapper.py --input1=$i --input2=${input2}.$splitno --outsuffix=.${splitno} $* > /mnt/galaxyData/tmp/palmapper.$suffix.$splitno.log" | ${qsub} -cwd -N "plm$suffix$splitno" -pe "*" 1-4 -j y -o /mnt/galaxyData/tmp/palmapper.$suffix.$splitno.log
    fi
done

count_running=`$qstat | grep $suffix | wc -l`
count_running_all=`$qstat | grep $suffix | wc -l`
while [ "$count_running" != 0 ];
do
    count_running=`$qstat | grep $suffix | wc -l`
    echo `date`: jobs still running: $count_running/$count_running_all
    sleep 10
done

samtoolspar=
junctionsoutput=
bamoutput=
samoutput=
logoutput=

while [ "$#" != "0" ]; do
    if [ "${1:0:14}" == "--bamsort=name" ];
    then
	samtoolpar='-n'
    fi
    if [ "${1:0:13}" == "--bam-output=" ];
    then
	bamoutput=${1:13:1000}
    fi
    if [ "${1:0:13}" == "--sam-output=" ];
    then
	samoutput=${1:13:1000}
    fi
    if [ "${1:0:19}" == "--junctions-output=" ];
    then
	junctionsoutput=${1:19:1000}
    fi
    if [ "${1:0:10}" == "--logfile=" ];
    then
	logoutput=${1:10:1000}
    fi
    shift
done


if [ "$bamoutput" != "None" ]; then
    echo merging bam file
    samtools merge $samtoolspar $bamoutput $bamoutput.????
    rm -f $bamoutput.????
fi

if [ "$samoutput" != "None" ]; then
    echo merging sam file

    cat /dev/null > $samoutput

    for i in ${input1}.????; do
	splitno=${i:(-4)} 
	cat ${samoutput}.$splitno >> $samoutput
	rm -f ${samoutput}.$splitno
    done
fi

if [ "$junctionsoutput" != "None" ]; then
    echo merging junction file

    cat /dev/null > $junctionsoutput

    for i in ${input1}.????; do
	splitno=${i:(-4)} 
	cat ${junctionsoutput}.$splitno >> $junctionsoutput
	rm -f ${junctionsoutput}.$splitno
    done
fi

if [ "$logoutput" != "None" ]; then

    for i in ${input1}.????; do
	splitno=${i:(-4)} 
	cat ${logoutput}.$splitno >> $logoutput
	cat /mnt/galaxyData/tmp/palmapper.$suffix.$splitno.log >> $logoutput
	rm -f ${logoutput}.$splitno /mnt/galaxyData/tmp/palmapper.$suffix.$splitno.log
    done
fi

rm -f ${input1}.????

echo done

date
