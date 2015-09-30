/*
 * PATH_WEIGHTS.CPP	
 *	    
 * The calling syntax is:
 *
 * [trans_weights, plif_weights] = path_weights(state_seq, obs_seq, score_plifs, num_states)
 *
 * Counts how often individual transitions/supporting points of feature
 * scoring functions are used by the given state sequence.
 *
 * state_seq    -- sequence of states for which weights are computed
 * obs_seq      -- sequence of observations, i.e. the feature matrix
 * score_plifs  -- a struct representation of feature scoring functions
 *                 (see also score_plif_struct.h / .cpp)
 * num_states   -- the number of states (in state_model)
 *
 * returns weights of transitions, i.e. counts of how often
 *   a certain transition has been used for decoding (trans_weights)
 *   and weights of supporting points of the feature scoring functions
 *   indicating how often certain scores are used (plif_weights)
 *
 * Compile using
 *   mex path_weights.cpp score_plif_struct.cpp 
 *
 * Written by Georg Zeller, MPI Tuebingen, Germany
 */

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "mex.h"
//#include <matrix.h>
#include "score_plif_struct.h"


void add_fraction(double* vec, const double val, const double* lim, const int L) {
  int idx = 0;
  for (int i=0; i<L; ++i) {
    if (val < lim[i])
      break;
    
    idx = idx + 1;
  }

  if (idx == 0) {
    ++vec[0];
  } else if (idx == L) {
      ++vec[L-1];
  } else {
    const double diff = lim[idx] - lim[idx-1];
    vec[idx]   +=  (val - lim[idx-1]) / diff;  
    vec[idx-1] += (lim[idx] - val) / diff;
  }
}

/*
 * index conversion from a pair of (row, column) indices (Matlab
 * style) into a single array index (C style)
 */
inline int conv_index2(const int row, const int col, const int NUM_ROWS) {
  return col * NUM_ROWS + row;
}

/*
 * index conversion from a triple of (d1, d2, d3) indices (Matlab
 * style for 3-dimensional matrices) into a single array index (C
 * style)
 */
inline int conv_index3(const int d1, const int d2, const int d3, const int DIM1, const int DIM2) {
  return d3 * DIM1*DIM2 + d2*DIM1 + d1;
}


/* actual mex function for computing path weights */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if(nrhs != 4) {
    mexErrMsgTxt("expected 4 input arguments:\n [trans_weights, plif_weights] = path_weights(state_seq, obs_seq, score_plifs, num_states)\n");
    return;
  }
  
  if(nlhs != 2) {
    mexErrMsgTxt("expected 2 output argument\n [trans_weights, plif_weights] = path_weights(state_seq, obs_seq, score_plifs, num_states)\n");
    return;
  }
  
//  mexPrintf("reading input arguments...\n");
  int arg = 0;

  const int m = mxGetM(prhs[arg]);
  const int LEN = mxGetN(prhs[arg]);
  assert(m == 1);
  double *state_seq = mxGetPr(prhs[arg]);
  ++arg;

  const int NUM_FEATS = mxGetM(prhs[arg]);
  const int n = mxGetN(prhs[arg]);
  assert(LEN == n);
  double *obs_seq = mxGetPr(prhs[arg]);
  ++arg;

  score_plif_struct *scr_ptr;
  int scr_M = 0; int scr_N = 0;
  scr_ptr = read_score_plif_struct(prhs[arg], scr_M, scr_N);
  //mexPrintf("  read score PLiF struct (%i x %i)\n", scr_M, scr_N);
  ++arg;

  const int NUM_STATES = (int) mxGetScalar(prhs[arg]);
  //mexPrintf("finished reading input data\n");

  // create 1st return argument: array of transition weights
  plhs[0] = mxCreateDoubleMatrix(NUM_STATES, NUM_STATES, mxREAL);
  double *trans_weights = mxGetPr(plhs[0]);
  // will already be intialized to 0 by mxCreateDoubleMatrix

  // compute transition weights
  for (int i=0; i<LEN-1; ++i) {
    // state indices are 1-based in Matlab, but here have to be 0-based
    const int s1 = state_seq[i] - 1;
    const int s2 = state_seq[i+1] - 1;
    const int p = conv_index2(s1, s2, NUM_STATES);
    ++trans_weights[p];
  }


  // create 2nd return argument: 3-dim. array of plif weights
  const int num_dim = 3;
  const int NUM_PLIF_NODES = scr_ptr[0].len;
  mwSize* dims;
  dims = (mwSize*) mxMalloc (num_dim * sizeof(mwSize));
  //int dims[3];
  dims[0] = NUM_FEATS;
  dims[1] = NUM_STATES;
  dims[2] = NUM_PLIF_NODES;
  plhs[1] = mxCreateNumericArray(num_dim,dims, mxDOUBLE_CLASS, mxREAL);
  double *plif_weights = mxGetPr(plhs[1]);
  // will already be intialized to 0 by mxCreateNumericArray

  // compute plif weights
  double vec[NUM_PLIF_NODES];
  int vec_idx[NUM_PLIF_NODES];
  for (int f=0; f<NUM_FEATS; ++f) { // for all features
    for (int p=0; p<LEN; ++p) {     // for all positions in the sequence
      const int s = state_seq[p] - 1;
      const score_plif_struct scr = scr_ptr[conv_index2(f, s, NUM_FEATS)];
      const double val = obs_seq[conv_index2(f, p, NUM_FEATS)];
      //mexPrintf("f = %i, p = %i, s = %i\n", f, p, s);
      for (int n=0; n<NUM_PLIF_NODES; ++n) {
        vec_idx[n] = conv_index3(f, s, n, NUM_FEATS, NUM_STATES);
        vec[n] = plif_weights[vec_idx[n]];
      }
      add_fraction(vec, val, scr.limits, NUM_PLIF_NODES);
      for (int n=0; n<NUM_PLIF_NODES; ++n) {
        plif_weights[vec_idx[n]] = vec[n];
        //mexPrintf("vec[%i] = %f (idx = %i)\n", n, vec[n], vec_idx[n]);
      }
    }
  }

  // clean up
  delete_score_plif_struct_matrix(scr_ptr, scr_M, scr_N);
}

