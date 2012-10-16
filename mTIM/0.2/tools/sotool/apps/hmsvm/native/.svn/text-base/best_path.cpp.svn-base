/*
 * BEST_PATH.CPP
 *     computes the state sequence maximizing a given scoring scheme 
 *     (precomputed "emission" scores score_matrix, transition scores A, 
 *     a distribution of start states as well as end states)
 *     using the Viterbi algorithm
 *
 * The calling syntax is:
 *
 *		[score state_seq] = best_path(p, q, A, score_matrix) 
 *
 * Compile using
 *   mex best_path.cpp
 *
 * Written by Gunnar Raetsch, Soeren Sonnenburg & Georg Zeller, MPI Tuebingen, Germany
 * Adpoted from the Shogun function best_path_trans_simple (see www.shogun-toolbox.org)
 */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "mex.h"
#include "score_plif_struct.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {		  

  if(nrhs != 4) {
    mexErrMsgTxt("expected 4 input arguments:\n [score state_seq] = best_path(p, q, A, score_matrix)\n");
    return;
  }
  
  if(nlhs != 2) {
    mexErrMsgTxt("expected 2 output argument2\n [score state_seq] = best_path(p, q, A, score_matrix)\n");
    return;
  }
  
//  fprintf(stderr, "reading input arguments...\n");
  
  // read input arguments
  int arg = 0;

  // distribution of allowed start states
  double *p;
  const int M_p = mxGetM(prhs[arg]);
  const int N_p = mxGetN(prhs[arg]);
  p = mxGetPr(prhs[arg]);
  ++arg;
  
  // distribution of allowed end states
  double *q;
  const int M_q = mxGetM(prhs[arg]);
  const int N_q = mxGetN(prhs[arg]);
  q = mxGetPr(prhs[arg]);
  ++arg;
  
  // square matrix of allowed transitions (and their scores) 
  double *A;
  const int M_A = mxGetM(prhs[arg]);
  const int N_A = mxGetN(prhs[arg]);
  A = mxGetPr(prhs[arg]);
  ++arg;

  // matrix of "emission" scores
  double *scr;
  const int M_scr = mxGetM(prhs[arg]);
  const int N_scr = mxGetN(prhs[arg]);
  scr = mxGetPr(prhs[arg]);
  ++arg;

  const int LEN = N_scr;        // length of the given example
  const int NUM_STATES = M_scr; // number of states
  //assert(M_p == NUM_STATES);
  //assert(M_q == NUM_STATES);
  //assert(M_A == NUM_STATES);
  //assert(N_A == NUM_STATES);

  // start dynamic programming
//  printf("filling dynamic programmming matrix...\n");
//  printf("NUM_STATES=%i, LEN=%i\n", NUM_STATES, LEN);
  const double INF = INFINITY;
  double **dpm = new double*[LEN];   // d.p. matrix
  if (dpm==NULL)
    mexErrMsgTxt("memory allocation failed");

  int **trb = new int*[LEN];         // traceback matrix
  if (trb==NULL)
    mexErrMsgTxt("memory allocation failed");
  
  for (int i = 0; i<LEN; ++i) {
    dpm[i] = new double[NUM_STATES];
    if (dpm[i]==NULL)
      mexErrMsgTxt("memory allocation failed");

    trb[i] = new int[NUM_STATES];
    if (trb[i]==NULL)
      mexErrMsgTxt("memory allocation failed");
  }

  for (int s=0; s<NUM_STATES; ++s) {
    if (p[s] > -INF) {
      // dpm[0][s] = scr[s][0];
      dpm[0][s] = scr[s];
      // printf("dpm[0][%i]=%f\n", s, dpm[0][s]);
    } else {
      dpm[0][s] = -INF;
      // printf("dpm[0][%i]=%f\n", s, dpm[0][s]);
    }
  }

  for (int p=1; p<LEN; ++p) {
    // is there at least one transition possible?
    int trans_possible = -1;

    for (int t=0; t<NUM_STATES; ++t) {
      dpm[p][t] = -INF;
      trb[p][t] = -1;
      // precomputed emission score e = scr[t][p]
      double e = scr[p*NUM_STATES+t];
      // printf("e(%i,%i)=%f\n", t, p, e);
      // find maximally scoring prefix
      for (int s=0; s<NUM_STATES; ++s) {
        // transition score A[s][t]
        double a = A[t*NUM_STATES+s];
        //printf("pos=%i:  from %i to %i a=%f\n", p, t, s, a);
        if (a>-INF) {
          // tmp_score = e + A[s][t] + dpm[p-1][s];
          double tmp_score = e + a + dpm[p-1][s];
          trans_possible = 1;
          if (tmp_score > dpm[p][t]) { 
            dpm[p][t] = tmp_score;
            trb[p][t] = s;
          }
        }
      }
    }
    // check if at least one transition is possible
    if (trans_possible == -1) {
      printf("pos=%i: All transitions from had value -inf\n", p);
      mexErrMsgTxt("Invalid transition scores!");
    }
    
  }
  
  // traceback
//  fprintf(stderr, "tracing back best path...\n");
  int *opt_path = new int[LEN];
  if (opt_path==NULL)
    mexErrMsgTxt("memory allocation failed");
  double opt_score = -INF;
  opt_path[LEN-1] = -1;
  for (int s=0; s<NUM_STATES; ++s) {
    if (q[s] > -INF && dpm[LEN-1][s] > opt_score) {
      opt_score = dpm[LEN-1][s];
      opt_path[LEN-1] = s;
    }
  }
  if (opt_path[LEN-1] == -1) {
    for (int s=0; s<NUM_STATES; ++s) {
      fprintf(stdout, "s=%i, dpm[LEN-1][s]=%f\n", s, dpm[LEN-1][s]);
    }
    mexErrMsgTxt("Error: no entry point for trace-back!");      
  }

  for (int p=LEN-1; p>0; --p) {
    if (trb[p][opt_path[p]] == -INF) {
	mexErrMsgTxt("Error: stuck in trace-back!");      
    }
    opt_path[p-1] = trb[p][opt_path[p]];
  }

  //assert(p[opt_path[1]] > -INF);
  //assert(q[opt_path[LEN-1]] > -INF);

  // prepare return values
//  fprintf(stderr, "writing return values...\n");

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  if (plhs[0]==NULL)
    mexErrMsgTxt("memory allocation failed");
  double *ret_score = mxGetPr(plhs[0]);
  *ret_score = opt_score;
  
  plhs[1] = mxCreateDoubleMatrix(1, LEN, mxREAL);
  if (plhs[1]==NULL)
    mexErrMsgTxt("memory allocation failed");
  double *ret_path = mxGetPr(plhs[1]);
  // convert 0-based state sequence to 1-based sequence for Matlab
  for (int p=0; p<LEN; ++p) {
//    printf("best_path[%i] = %i\n", p, opt_path[p]);
    ret_path[p] = ((double) opt_path[p]) + 1.0;
  }

  // clean-up
  for (int i = 0; i<LEN; ++i) {
    delete[] dpm[i];
    delete[] trb[i];
  }
  delete[] dpm;
  delete[] trb;
  delete[] opt_path;
}
