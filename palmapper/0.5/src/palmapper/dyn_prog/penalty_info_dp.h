// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#ifndef __PENALTY_INFO_DP_H__
#define __PENALTY_INFO_DP_H__

#include <assert.h>

#include "common_dp.h"
#include "Mathmatics_dp.h"

enum ETransformType
{
	T_LINEAR,
	T_LOG,
	T_LOG_PLUS1,
	T_LOG_PLUS3,
	T_LINEAR_PLUS3
}  ;

struct penalty_struct
{
	INT len ;
	REAL *limits ;
	REAL *penalties ;
	INT max_len ;
	INT min_len ;
	REAL *cache ;
	enum ETransformType transform ;
	INT id ;
	struct penalty_struct *next_pen ;
	char * name ;
	INT use_svm ;
} ;

void init_penalty_struct(struct penalty_struct &PEN) ;
void delete_penalty_struct_palma(struct penalty_struct &PEN) ;
void init_penalty_struct_cache(struct penalty_struct &PEN) ;
//void delete_penalty_struct_array(struct penalty_struct *PEN, INT len) ;
void copy_penalty_struct(struct penalty_struct *old, struct penalty_struct *newp);

//#ifdef HAVE_MATLAB
//struct penalty_struct * read_penalty_struct_from_cell(const mxArray * mx_penalty_info, mwSize &P) ;
//#endif

REAL lookup_penalty(const struct penalty_struct *PEN, INT p_value, 
					REAL* svm_values, bool follow_next=true) ;

inline REAL lookup_penalty_dummy(const struct penalty_struct *PEN, INT p_value, 
						  REAL* svm_values, bool follow_next)
{
	if (PEN==NULL)
		return 0 ;
	assert(!PEN->use_svm) ;
		
	if ((p_value<PEN->min_len) || (p_value>PEN->max_len))
		return -ALMOST_INFINITY ;

	return 0;
}


enum mode { NORMAL, USE_QUALITY_SCORES};

#endif
