include ../../bin/mTIM_config.sh

MKOCTFILE=${OCTAVE_MKOCT}
MEX=${MATLAB_MEX_PATH}

all:	convert_states2labels.mex convert_states2labels.mexa64
octave:	convert_states2labels.mex 
matlab:	convert_states2labels.mexa64

convert_states2labels.mexa64:	convert_states2labels.cpp
	rm -f *.o
	${MEX} -g -O convert_states2labels.cpp 

convert_states2labels.mex:	convert_states2labels.cpp 
	rm -f *.o
	${MKOCTFILE} -g --mex convert_states2labels.cpp 

clean:	
	rm -f *.o *.mexa64 *.mex 
