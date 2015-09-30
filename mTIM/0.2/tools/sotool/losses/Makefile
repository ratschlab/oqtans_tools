include ./../../../bin/mTIM_config.sh

MEX=${MATLAB_MEX_PATH}
MKOCTFILE=${OCTAVE_MKOCT}

all: loss_hinge_native.mex loss_hinge_native.mexa64 loss_sqrhinge_native.mex loss_sqrhinge_native.mexa64 loss_logistic_native.mex loss_logistic_native.mexa64 
octave: loss_hinge_native.mex loss_sqrhinge_native.mex loss_logistic_native.mex
matlab: loss_hinge_native.mexa64 loss_sqrhinge_native.mexa64 loss_logistic_native.mexa64

loss_hinge_native.mex: loss_hinge_native.cpp 
	rm -f *.o
	${MKOCTFILE} -g --mex loss_hinge_native.cpp 

loss_hinge_native.mexa64: loss_hinge_native.cpp 
	rm -f *.o
	${MEX} -g -O loss_hinge_native.cpp  

loss_sqrhinge_native.mex: loss_sqrhinge_native.cpp
	rm -f *.o
	${MKOCTFILE} -g --mex loss_sqrhinge_native.cpp 

loss_sqrhinge_native.mexa64: loss_sqrhinge_native.cpp
	rm -f *.o
	${MEX} -g -O loss_sqrhinge_native.cpp

loss_logistic_native.mex: loss_logistic_native.cpp
	rm -f *.o
	${MKOCTFILE} -g --mex loss_logistic_native.cpp 

loss_logistic_native.mexa64: loss_logistic_native.cpp
	rm -f *.o
	${MEX} -g -O loss_logistic_native.cpp

clean:	
	rm -f *.o *.mexa64 *.mex
