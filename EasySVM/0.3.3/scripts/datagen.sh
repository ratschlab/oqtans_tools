#!/bin/bash
export PYTHONPATH=/mnt/oqtansTools/oqtans_dep/lib/python2.7/site-packages/:/mnt/oqtansTools/oqtans_dep/lib/python2.7/dist-packages/:$PYTHONPATH
export LD_LIBRARY_PATH=/mnt/oqtansTools/oqtans_dep/lib/:$LD_LIBRARY_PATH
/usr/bin/python -W ignore::UserWarning /mnt/oqtansTools/oqtans/EasySVM/0.3.3/scripts/datagen.py $@
