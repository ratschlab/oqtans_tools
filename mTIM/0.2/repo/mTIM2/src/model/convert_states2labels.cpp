/*
 * CONVERT_STATES2LABELS.CPP	
 *	    converts a state sequence into a label sequence
 *
 * The calling syntax is:
 *
 *		label_seq = convert_states2labels(state_seq, state_labeling) 
 *
 * Compile using
 *   mex convert_states2labels.cpp
 *
 * Written by Georg Zeller, MPI Tuebingen, Germany
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {		  

  if(nrhs != 2) {
    mexErrMsgTxt("expected 2 input arguments:\n label_seq = convert_states2labels(state_seq, state_labeling)\n");
    return;
  }
  
  if(nlhs != 1) {
    mexErrMsgTxt("expected 1 output argument\n label_seq = convert_states2labels(state_seq, state_labeling)\n");
    return;
  }
  
// printf("reading input arguments...\n");
  int arg = 0;

  const int m = mxGetM(prhs[arg]);
  const int L = mxGetN(prhs[arg]);
  assert(m == 1);
  double *state_seq = mxGetPr(prhs[arg]);
  ++arg;

  const int S = mxGetM(prhs[arg]);
  const int n = mxGetN(prhs[arg]);
  // assert(n == 1);
  double *state_labeling = mxGetPr(prhs[arg]);
	
// printf("finished reading input data\n");

  //fprintf(stdout, "S=%i, L=%i, n=%i\n", S, L, n);
  // S is the number of states
  // L is the length of the sequence
  // create label sequence which will be returned
  plhs[0] = mxCreateDoubleMatrix(1, L, mxREAL);
  double *label_seq = mxGetPr(plhs[0]);

  for (int pos=0; pos<L; ++pos) {
    const int state = ((int) state_seq[pos]) - 1;
    label_seq[pos] = state_labeling[state];
  }
}
