/*
 * COMPUTE_SCORE_MATRIX.CPP	
 *	    computes a matrix for score cumulation during decoding 
 *          given a training example X, and feature scoring functions score_fcts
 *
 * The calling syntax is:
 *
 *		[scr_matrix] = compute_score_matrix(X,score_fcts) 
 *
 * Compile using
 *   mex compute_score_matrix.cpp score_plif_struct.cpp 
 *
 * Written by Gunnar Raetsch & Georg Zeller, MPI Tuebingen, Germany
 */

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//#include "matrix.h"
#include "mex.h"
#include "score_plif_struct.h"

/*
 * index conversion from a pair of (row, column) indices (Matlab
 * style) into a single array index (C style)
 */
inline int conv_index2(const int row, const int col, const int NUM_ROWS) {
  return col * NUM_ROWS + row;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {		  

  if(nrhs < 2 || nrhs==3 || nrhs > 4) {
    mexErrMsgTxt("expected either 2 or 4 input arguments:\n [scr_matrix] = compute_score_matrix(X,score_fcts[,min_feat_values,max_feat_values])\n");
    return;
  }
  
  if(nlhs != 1) {
    mexErrMsgTxt("expected either 2 or 4 input arguments:\n [scr_matrix] = compute_score_matrix(X,score_fcts[,min_feat_values,max_feat_values])\n");
    return;
  }
  
  // read input arguments
  int arg = 0;

  double *X_ptr;
  const int X_M = mxGetM(prhs[arg]);
  const int X_N = mxGetN(prhs[arg]);
  X_ptr = mxGetPr(prhs[arg]);
  ++arg;
  
  score_plif_struct *scr_ptr;
  int scr_M = 0; int scr_N = 0;
  scr_ptr = read_score_plif_struct(prhs[arg], scr_M, scr_N);
  ++arg;

  const double INF = mxGetInf();
  //const double INF = INFINITY;
  const int L = X_N;                      // length of the given example
  const int n_states = scr_N;             // number of states
  const int n_feats = scr_M;              // number of features

  double *min_feat_val;
  double *max_feat_val;
  if (nrhs == 4) {
    min_feat_val = mxGetPr(prhs[arg]);
    ++arg;
    
    max_feat_val = mxGetPr(prhs[arg]);
    ++arg;
  } else {
    // intialize min_feat_val and max_feat_val to -inf and +inf, respectively
    mxArray *mn = mxCreateDoubleMatrix(n_states, n_feats, mxREAL);
    min_feat_val = mxGetPr(mn);
    mxArray *mx = mxCreateDoubleMatrix(n_states, n_feats, mxREAL);
    max_feat_val = mxGetPr(mx);
    for (int f=0; f<n_feats; ++f) {       // for all features
      for (int s=0; s<n_states; ++s) {    // for all states
	    min_feat_val[f*n_states+s] = -INF;
	    max_feat_val[f*n_states+s] =  INF;
      }
    }
  }

  // compute scores along the possible paths for X
  plhs[0] = mxCreateDoubleMatrix(n_states, L, mxREAL);
  double *pp = mxGetPr(plhs[0]);
  // will already be initialized to 0.0 by mxCreateDoubleMatrix

  // scores for real-valued features
  for (int pos=0; pos<L; ++pos) {         // for all positions in given example
    for (int f=0; f<n_feats; ++f) {       // for all features
      for (int s=0; s<n_states; ++s) {    // for all states
	    const int scr_idx = s*n_feats+f;
	    const int pp_idx = pos*n_states+s;
	    const int X_idx = pos*X_M+f;

        /*
        if (pos == 0) {
          printf("f=%i, s=%i, min_feat_val[s,f]=%f, max_feat_val[s,f]=%f\n", f, s, min_feat_val[f*n_states+s], max_feat_val[f*n_states+s]);
        }
        */

        if (X_ptr[X_idx] >= min_feat_val[f*n_states+s] && X_ptr[X_idx] <= max_feat_val[f*n_states+s]) {
          pp[pp_idx] += lookup_score_plif(&scr_ptr[scr_idx], X_ptr[X_idx]);
          //	  if (pos == 0)
          //	    printf("  permitted\n");
        } else {
          pp[pp_idx] = -INF;
        }
      }
    }
  }
  	
  // clean up
  delete_score_plif_struct_matrix(scr_ptr, scr_M, scr_N);
}
