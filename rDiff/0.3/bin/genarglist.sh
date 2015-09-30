#/bin/bash

list=
until [ -z $1 ] ; do
	if [ $# != 1 ];
	then
		list="${list}$1:"
	else
		list="${list}$1"
	fi
	shift
done
echo $list

