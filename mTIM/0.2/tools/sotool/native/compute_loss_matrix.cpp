/*
 * COMPUTE_LOSS_MATRIX.CPP	
 *	    computes a loss matrix for decoding the max-margin violator 
 *
 * The calling syntax is:
 *
 *		[loss_matrix] = compute_loss_matrix(loss,true_states) 
 *
 * Compile using
 *   mex compute_loss_matrix.cpp
 *
 * Written by Georg Zeller, MPI Tuebingen, Germany
 */

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {		  

  if(nrhs != 2) {
    mexErrMsgTxt("expected 2 input arguments:\n [loss_matrix] = compute_loss_matrix(loss,true_states)\n");
    return;
  }
  
  if(nlhs != 1) {
    mexErrMsgTxt("expected 1 output argument\n [loss_matrix] = compute_loss_matrix(loss,true_states)\n");
    return;
  }
  
// printf("reading input arguments...\n");
  int arg = 0;

  const int M = mxGetM(prhs[arg]);
  const int n = mxGetN(prhs[arg]);
  assert(M == n);
  double *loss = mxGetPr(prhs[arg]);
  ++arg;
	
  const int m = mxGetM(prhs[arg]);
  const int L = mxGetN(prhs[arg]);
  assert(m == 1);
  double *state_seq = mxGetPr(prhs[arg]);
//  printf("finished reading input data\n");

// printf("M=%i, L=%i\n", M, L);
  // M is the number of states
  // L is the length of the sequence
  // create loss matrix which will be returned
  plhs[0] = mxCreateDoubleMatrix(M, L, mxREAL);
  double *lm = mxGetPr(plhs[0]);

  for (int i=0; i<L*M; ++i)
    lm[i] = 0.0;
  // fill loss matrix
  for (int pos=0; pos<L; ++pos) { // for all positions in given example
    for (int s=0; s<M; ++s) {     // for all states
	const int lm_idx = pos*M+s;
	const int true_state = ((int) state_seq[pos]) - 1;
	lm[lm_idx] = loss[true_state*M+s];
    }
  }
}
