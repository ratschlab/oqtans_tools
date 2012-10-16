#ifndef __SCORE_PLIF_STRUCT_H__
#define __SCORE_PLIF_STRUCT_H__

#include <mex.h>

/*
 * Defines a struct representation of feature scoring functions
 *
 * Written by Gunnar Raetsch & Georg Zeller, MPI Tuebingen, Germany 
 */

struct score_plif_struct
{
	int len;
	double *limits;
	double *scores;
	double *weights;
	int feat_idx;
};

void init_score_plif_struct(struct score_plif_struct &SCR);
void delete_score_plif_struct_matrix(struct score_plif_struct *SCR, int &M, int &N);

struct score_plif_struct* read_score_plif_struct(const mxArray *mx_score_plif_matrix, int &M, int &N);
mxArray* write_score_plif_struct(const struct score_plif_struct *SCR, int &M, int &N);
#endif

double lookup_score_plif(const struct score_plif_struct *SCR, double value);
void update_score_plif_weights(struct score_plif_struct *SCR, const double value);
