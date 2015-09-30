#!/bin/bash
#
# oqtans sub module compile program 
#
set -e 

source ./../../../oqtans_config.sh

cd $OQTANS_PATH/mTIM/0.2/src/utils
if [ "$1" == "" -o "$1" == "all" ];
then
    make octave
    cd $OQTANS_PATH/mTIM/0.2/tools/sotool/losses
    make octave 
    cd $OQTANS_PATH/mTIM/0.2/tools/sotool/native
    make octave
    cd $OQTANS_PATH/mTIM/0.2/src/model 
    make octave
fi

if [ "$1" == "clean" ];
then
    make clean
    cd $OQTANS_PATH/mTIM/0.2/tools/sotool/losses
    make clean
    cd $OQTANS_PATH/mTIM/0.2/tools/sotool/native
    make clean
    cd $OQTANS_PATH/mTIM/0.2/src/model 
    make clean
fi
