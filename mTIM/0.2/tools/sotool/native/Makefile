include ./../../../bin/mTIM_config.sh

MEX=${MATLAB_MEX_PATH}
MKOCTFILE=${OCTAVE_MKOCT}

all: compute_score_matrix.mex compute_score_matrix.mexa64 best_path.mex best_path.mexa64 compute_loss_matrix.mex compute_loss_matrix.mexa64 path_weights.mex path_weights.mexa64 
octave: compute_score_matrix.mex best_path.mex compute_loss_matrix.mex path_weights.mex
matlab: compute_score_matrix.mexa64 best_path.mexa64 compute_loss_matrix.mexa64 path_weights.mexa64

compute_score_matrix.mex: compute_score_matrix.cpp score_plif_struct.cpp
	rm -f *.o
	${MKOCTFILE} -g --mex compute_score_matrix.cpp score_plif_struct.cpp 

compute_score_matrix.mexa64: compute_score_matrix.cpp score_plif_struct.cpp
	rm -f *.o
	${MEX} -g -O compute_score_matrix.cpp score_plif_struct.cpp

best_path.mex: best_path.cpp 
	rm -f *.o
	${MKOCTFILE} -g --mex best_path.cpp 

best_path.mexa64: best_path.cpp
	rm -f *.o
	${MEX} -g -O best_path.cpp 

compute_loss_matrix.mex: compute_loss_matrix.cpp 
	rm -f *.o
	${MKOCTFILE} -g --mex compute_loss_matrix.cpp 

compute_loss_matrix.mexa64: compute_loss_matrix.cpp 
	rm -f *.o
	${MEX} -g -O compute_loss_matrix.cpp

path_weights.mex: path_weights.cpp score_plif_struct.cpp  
	rm -f *.o
	${MKOCTFILE} -g --mex path_weights.cpp score_plif_struct.cpp 

path_weights.mexa64: path_weights.cpp score_plif_struct.cpp  
	rm -f *.o
	${MEX} -g -O path_weights.cpp score_plif_struct.cpp

clean:	
	rm -f *.o *.mexa64 *.mex
