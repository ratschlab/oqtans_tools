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

set -e 

PROG=`basename $0`
DIR=`dirname $0`

ORGANISM=${1}
GENOME_SEQ=${2}
GFF3_INPUT=${3}
SAM_INPUT=${4}
INFO=${5}
PROFILES_FN=${6}

. ${DIR}/../bin/rquant_config.sh
${DIR}/../bin/rquant_gendata ${ORGANISM} ${GENOME_SEQ} ${GFF3_INPUT} ${SAM_INPUT} ${INFO} ${PROFILES_FN}
