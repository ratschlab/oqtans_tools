#/bin/bash
# Copyright (C) 2010-2012 Max Planck Society

until [ -z $1 ] ; do
	if [ $# != 1 ];
	then
		echo -n "'$1', "
	else
		echo -n "'$1'"
	fi
	shift
done
