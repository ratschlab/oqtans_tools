#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "mex.h"

#include "score_plif_struct.h"

/*
 * Defines a struct representation of feature scoring functions
 *
 * Written by Gunnar Raetsch & Georg Zeller, MPI Tuebingen, Germany 
 */

inline bool my_isnan(double x) {
  return x != x;
}

void init_score_plif_struct(struct score_plif_struct &SCR) {
  SCR.limits = NULL;
  SCR.scores = NULL;
  SCR.weights = NULL;
  SCR.feat_idx = -1;
}

void delete_score_plif_struct_matrix(struct score_plif_struct *SCR, int &M, int &N) {
  for (int i=0; i<M*N; ++i) {
    if (SCR[i].feat_idx!=-1) {
      delete[] SCR[i].limits;
      delete[] SCR[i].scores;
      delete[] SCR[i].weights;
    }
  }
  delete[] SCR;
}


struct score_plif_struct* read_score_plif_struct(const mxArray *mx_score_plif_matrix, int &M, int &N) {
  int n_fields = mxGetNumberOfFields(mx_score_plif_matrix);
  int n_elems  = mxGetNumberOfElements(mx_score_plif_matrix);
//  fprintf(stderr, "n_fields %i, n_elems %i\n", n_fields, n_elems);

  M = mxGetM(mx_score_plif_matrix);
  N = mxGetN(mx_score_plif_matrix);
	
//  fprintf(stderr, "Initializing score_plif struct matrix of dimensions %i x %i\n", M, N) ;
  struct score_plif_struct *SCR = new struct score_plif_struct[M*N];
  for (int i=0; i<M; ++i) {
    for(int j=0; j<N; ++j) {
      init_score_plif_struct(SCR[i*N+j]);
//	fprintf(stderr, "initialized field (%i,%i)\n", i, j);
    }
  }
  
//  fprintf(stderr, "Reading limits...\n") ;
  int lim_f_num = mxGetFieldNumber(mx_score_plif_matrix, "limits");
  int scr_f_num = mxGetFieldNumber(mx_score_plif_matrix, "scores");
  int dim_f_num = mxGetFieldNumber(mx_score_plif_matrix, "dim");
//  fprintf(stderr, "%i  %i  %i", lim_f_num, scr_f_num, dim_f_num) ;
	
  for(int j=0; j<n_elems; ++j) {
    const mxArray *lim = mxGetFieldByNumber(mx_score_plif_matrix, j, lim_f_num);
    const mxArray *scr = mxGetFieldByNumber(mx_score_plif_matrix, j, scr_f_num);
    const mxArray *dim = mxGetFieldByNumber(mx_score_plif_matrix, j, dim_f_num);
//    fprintf(stderr, "%s  %s  %s", mxGetClassName(lim), mxGetClassName(scr),  mxGetClassName(dim)) ;

    if (lim==NULL || !mxIsDouble(lim)) {
      fprintf(stderr, "limits are expected to be doubles\n");
      delete_score_plif_struct_matrix(SCR,M,N);
      return NULL;
    }		
    
    if (scr==NULL || !mxIsDouble(scr)) {
      fprintf(stderr, "scores are expected to be doubles\n");
      delete_score_plif_struct_matrix(SCR,M,N);
      return NULL;
    }		
			 
    if (dim==NULL || !mxIsNumeric(dim)) {
      fprintf(stderr, "dim is expected to be a single (integer) value\n");
      delete_score_plif_struct_matrix(SCR,M,N);
      return NULL;
    }		
	
    const double *lim_ptr = mxGetPr(lim);
    const double *scr_ptr = mxGetPr(scr);
    const int feat = (int) *mxGetPr(dim);
    const int len = mxGetN(lim);
  			 
    SCR[j].len = len ;
    SCR[j].feat_idx = feat;
    SCR[j].limits = new double[len];
    SCR[j].scores = new double[len];
    for (int i=0; i<len; ++i) {
      SCR[j].limits[i]=lim_ptr[i];
      SCR[j].scores[i]=scr_ptr[i];
    }
  }
  return SCR;
}

double lookup_score_plif(const struct score_plif_struct *SCR, double value) {	
  if (SCR==NULL)
    return 0;
  
  int idx = 0;
  double ret;
  for (int i=0; i<SCR->len; i++)
    if (SCR->limits[i]<=value)
      idx++;
	
//  fprintf(stderr, "value:%f idx:%i limits[idx-1]:%f limits[idx]:%f\n", value, idx, SCR->limits[idx-1], SCR->limits[idx]) ;
	
  if (idx==0)
    ret=SCR->scores[0];
  else if (idx==SCR->len)
    ret=SCR->scores[SCR->len-1];
  else {
    ret = (SCR->scores[idx]*(value-SCR->limits[idx-1]) + SCR->scores[idx-1]*
	   (SCR->limits[idx]-value)) / (SCR->limits[idx]-SCR->limits[idx-1]);  
  }
  return ret;
}

void update_score_plif_weights(struct score_plif_struct *SCR, const double value) {
  if (SCR==NULL)
    return;
  
  int idx = -1;
  for (int i=0; i<SCR->len; ++i) {
    if (SCR->limits[i]<=value)
      ++idx;
    else
      break;
  }
//  fprintf(stderr, " idx %i limit[idx] %f limit[idx+1] %f value %f\n", idx, SCR->limits[idx], SCR->limits[idx+1], value);

  if (idx==-1)
    SCR->weights[0] += 1.0;
  else if (idx==SCR->len-1)
    SCR->weights[SCR->len-1] += 1.0;
  else {
    SCR->weights[idx+1] += (value - SCR->limits[idx]) / (SCR->limits[idx+1] - SCR->limits[idx]);
    SCR->weights[idx]   += (SCR->limits[idx+1] - value) / (SCR->limits[idx+1] - SCR->limits[idx]);
    if(my_isnan(SCR->weights[idx]) || my_isnan(SCR->weights[idx+1]))
      fprintf(stderr, " idx %i limit[idx] %f limit[idx+1] %f value %f\n", idx, SCR->limits[idx], SCR->limits[idx+1], value);
  }
}


mxArray* write_score_plif_struct(const struct score_plif_struct *SCR,  int &M, int &N) {
//  fprintf(stderr, "M: %i N: %i field_names ", M, N) ;
  const int nfields = 4;
  const char *field_names[] = {"limits", "scores", "dim", "w"};
  
  mxArray *ret = mxCreateStructMatrix(M, N, nfields, field_names);	
  int n_elems = mxGetNumberOfElements(ret);
  for(int j=0; j<n_elems; ++j) {
    int len = SCR[j].len;
    const int dims[] = {1,len};
    
    mxArray *lim = mxCreateNumericArray(2, (const mwSize*) dims, mxDOUBLE_CLASS, mxREAL);
    double *lim_ptr = mxGetPr(lim);
    memcpy(lim_ptr, SCR[j].limits, len*mxGetElementSize(lim));
    mxSetField(ret, j, field_names[0], lim);
    
    mxArray *scr = mxCreateNumericArray(2,(const mwSize*) dims, mxDOUBLE_CLASS, mxREAL);
    double *scr_ptr = mxGetPr(scr);
    memcpy(scr_ptr, SCR[j].scores, len*mxGetElementSize(scr));
    mxSetField(ret, j, field_names[1], scr);
    
    mxArray *dim = mxCreateDoubleScalar(SCR[j].feat_idx);
    mxSetField(ret, j, field_names[2], dim);
    
    mxArray *wght = mxCreateNumericArray(2, (const mwSize*) dims, mxDOUBLE_CLASS, mxREAL);
    double *w_ptr = mxGetPr(wght);
    memcpy(w_ptr, SCR[j].weights, len*mxGetElementSize(wght));
    mxSetField(ret, j, field_names[3], wght);
  }
  return ret;
}
