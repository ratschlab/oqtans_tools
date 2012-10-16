#!/bin/bash

#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# Written (W) 2009-2010 Regina Bohnert, Gunnar Raetsch
# Copyright (C) 2009-2010 Max Planck Society
#

# rQuant wrapper script to start the interpreter with the correct list of arguments

set -e

PROG=`basename $0`
DIR=`dirname $0`

exec ${DIR}/start_interpreter.sh ${PROG} "`${DIR}/genarglist.sh $@`"
