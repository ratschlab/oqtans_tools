MATLAB_DIR = /usr/local

MEX = $(MATLAB_DIR)/matlab2010a/bin/mex

MATLAB_INCL = $(MATLAB_DIR)/matlab2010a/extern/include

compute_score_matrix: 
	${MEX} compute_score_matrix.cpp score_plif_struct.cpp

best_path:
	${MEX} best_path.cpp

compute_loss_matrix: 
	${MEX} compute_loss_matrix.cpp

path_weights: 
	${MEX} path_weights.cpp score_plif_struct.cpp



all:	compute_score_matrix best_path compute_loss_matrix path_weights

clean:	
	rm -f *.o *.mexa64 *.mex
