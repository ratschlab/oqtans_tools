#!/bin/bash
# deseq-hts wrapper script to start the interpreter with the correct list of arguments
# Copyright (C) 2010-2012 Max Planck Society
set -e
PROG=`basename $0`
DIR=`dirname $0`
exec ${DIR}/start_interpreter.sh ${PROG} "`${DIR}/genarglist.sh $@`"
