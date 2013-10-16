/*
 * LOSS_HINGE_NATIVE.CPP	
 *
 * The calling syntax is:
 *
 *		[loss, idxs, params] = loss_hinge_native(t,delta_y_ybar,delta_psis,delta_psis_idxs) 
 *
 * Compile using
 *   mex loss_hinge_native.cpp
 *
 * Written by Nico Goernitz, TU Berlin, MPI Tuebingen, Germany
 */

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {		  

  if(nrhs < 5 || nrhs > 5) {
    mexErrMsgTxt("expected either 5 input arguments.\n");
		return;
  }
  
  if(nlhs != 3) {
    mexErrMsgTxt("expected 3 output arguments.\n");
    return;
  }
  
  // read input arguments
  int arg = 0;
  
  // amount of training examples == number of slacks
	const mxArray *params = prhs[arg];
	const int num_exm = (int) *mxGetPr(mxGetField(params,0, "num_examples"));
//	printf("Number of examples (new): %i\n",num_exm);
  ++arg;
	
	// this is the current solution
	const int M3 = mxGetM(prhs[arg]);
  const int D = mxGetN(prhs[arg]);
/*
  assert(M3 == 1);
*/
	double *w = mxGetPr(prhs[arg]);
	++arg;
	// get the delta_y_ybar [1 x t]-array
	const int M1 = mxGetM(prhs[arg]);
  const int N1 = mxGetN(prhs[arg]);
  assert(M1 == 1);
	double *delta_y_ybar = mxGetPr(prhs[arg]);
	++arg;

	const int DIM = mxGetM(prhs[arg]);
  const int N = mxGetN(prhs[arg]);
	double *delta_psis = mxGetPr(prhs[arg]);
	++arg;

	// get the delta_psis_idxs [1 x n]-array
	const int M2 = mxGetM(prhs[arg]);
  const int N2 = mxGetN(prhs[arg]);
  assert(M2 == 1);
	double *delta_psis_idxs = mxGetPr(prhs[arg]);
	++arg;

	 

  // these are the output values (already initialized)
	plhs[0] = mxCreateDoubleMatrix(1, num_exm, mxREAL);
  double *losses = mxGetPr(plhs[0]);

 	plhs[1] = mxCreateDoubleMatrix(1, num_exm, mxREAL);
  double *idxs = mxGetPr(plhs[1]);

	  // copy the parameter structure without making changes
	// to the elements
	plhs[2] = mxDuplicateArray(params);



	// for all dpsis
	for (int n=0; n<N; n++) {
		double loss = delta_y_ybar[n];
		const int ind = delta_psis_idxs[n];

		// for all dimensions in dpsis
		for (int d=0; d<DIM; d++) {
			loss -= w[d]*delta_psis[DIM*n+d];
		}
		
		// hinge 
		if (loss<0.0) {
			loss = 0.0;
		}

		//printf("t=%i n=%i loss=%f idxs=%i\n",num_exm,n,loss,ind);

		// store loss 
		if (loss>=losses[ind-1]) {
			losses[ind-1] = loss;
			idxs[ind-1] = n+1;
		}

	}

}
