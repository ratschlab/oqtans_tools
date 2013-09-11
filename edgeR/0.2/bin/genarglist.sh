#/bin/bash
##
# Copyright (C) 2009-2013 Max Planck Society & Memorial Sloan-Kettering Cancer Center
##

until [ -z $1 ] ; do
	if [ $# != 1 ];
	then
		echo -n "'$1', "
	else
		echo -n "'$1'"
	fi
	shift
done
