#!/bin/bash
until [ -z $1 ] ; do
	if [ $# != 1 ];
	then
		echo -n "'$1', "
	else
		echo -n "'$1'"
	fi
	shift
done
