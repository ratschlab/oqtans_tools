#/bin/bash
# rDiff wrapper script to start the interpreter with the correct list of arguments

set -e

PROG=`basename $0`
DIR=`dirname $0`

exec ${DIR}/start_interpreter.sh ${PROG} "`${DIR}/genarglist.sh $@`"
