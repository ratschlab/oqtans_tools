// Authors: Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean
// Copyright (C) 2006-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include <assert.h>
#include "config_dp.h"
//#include "features/CharFeatures.h"
//#include "features/StringFeatures.h"

#include <stdio.h>
#include <string.h>

#include "io_dp.h"


#include "penalty_info_dp.h"
#include "fill_matrix.h"

void init_penalty_struct(struct penalty_struct &PEN)
{
	PEN.limits=NULL ;
	PEN.penalties=NULL ;
	PEN.id=-1 ;
	PEN.next_pen=NULL ;
	PEN.transform = T_LINEAR ;
	PEN.name = NULL ;
	PEN.max_len=0 ;
	PEN.min_len=0 ;
	PEN.cache=NULL ;
	PEN.use_svm=0 ;
}

void init_penalty_struct_cache(struct penalty_struct &PEN)
{
	if (PEN.cache || PEN.use_svm)
		return ;
		
	REAL* cache=new REAL[PEN.max_len+1] ;
	if (cache)
	{
		for (INT i=0; i<=PEN.max_len; i++)
			cache[i] = lookup_penalty(&PEN, i, 0, false) ;
		PEN.cache = cache ;
	}
}

void delete_penalty_struct_palma(struct penalty_struct &PEN)
{
	if (PEN.id!=-1)
	{
		delete[] PEN.limits ;
		delete[] PEN.penalties ;
		delete[] PEN.name ;
		delete[] PEN.cache ;
	}
}

void delete_penalty_struct_array(struct penalty_struct *PEN, INT len)
{
	for (int i=0; i<len; i++)
		delete_penalty_struct_palma(PEN[i]) ;
	delete[] PEN ;
}

void copy_penalty_struct(struct penalty_struct *old, struct penalty_struct *newp) {
   newp->len       = old->len;
   newp->limits    = old->limits;
   newp->penalties = old->penalties;
   newp->max_len   = old->max_len;
   newp->min_len   = old->min_len;
   newp->cache     = old->cache;
   newp->transform = old->transform;
   newp->id        = old->id;
   newp->next_pen  = old->next_pen;
   newp->name      = old->name;
   newp->use_svm   = old->use_svm;
}

REAL lookup_penalty_svm(const struct penalty_struct *PEN, INT p_value, REAL *d_values)
{	
	if (PEN==NULL)
		return 0 ;
	assert(PEN->use_svm>0) ;
	REAL d_value=d_values[PEN->use_svm-1] ;
	//fprintf(stderr,"transform=%i, d_value=%1.2f\n", (INT)PEN->transform, d_value) ;
	
	switch (PEN->transform)
	{
	case T_LINEAR:
		break ;
	case T_LOG:
		d_value = log(d_value) ;
		break ;
	case T_LOG_PLUS1:
		d_value = log(d_value+1) ;
		break ;
	case T_LOG_PLUS3:
		d_value = log(d_value+3) ;
		break ;
	case T_LINEAR_PLUS3:
		d_value = d_value+3 ;
		break ;
	default:
		CIO_DP::message(M_ERROR_DP, "unknown transform\n") ;
		break ;
	}
	
	INT idx = 0 ;
	REAL ret ;
	for (INT i=0; i<PEN->len; i++)
		if (PEN->limits[i]<=d_value)
			idx++ ;
	
	if (idx==0)
		ret=PEN->penalties[0] ;
	else if (idx==PEN->len)
		ret=PEN->penalties[PEN->len-1] ;
	else
	{
		ret = (PEN->penalties[idx]*(d_value-PEN->limits[idx-1]) + PEN->penalties[idx-1]*
			   (PEN->limits[idx]-d_value)) / (PEN->limits[idx]-PEN->limits[idx-1]) ;  
	}
	
	//fprintf(stderr,"ret=%1.2f\n", ret) ;

	if (PEN->next_pen)
		ret+=lookup_penalty(PEN->next_pen, p_value, d_values);
	
	//fprintf(stderr,"ret=%1.2f\n", ret) ;
	
	return ret ;
}

REAL lookup_penalty(const struct penalty_struct *PEN, INT p_value, 
					REAL* svm_values, bool follow_next)
{	
	if (PEN==NULL)
		return 0 ;
	if (PEN->use_svm)
		return lookup_penalty_svm(PEN, p_value, svm_values) ;
		
	if ((p_value<PEN->min_len) || (p_value>PEN->max_len))
		return -ALMOST_INFINITY;
	
	if (PEN->cache!=NULL && (p_value>=0) && (p_value<=PEN->max_len))
	{
		REAL ret=PEN->cache[p_value] ;
		if (PEN->next_pen && follow_next)
			ret+=lookup_penalty(PEN->next_pen, p_value, svm_values);
		return ret ;
	}
	
	REAL d_value = (REAL) p_value ;
	switch (PEN->transform)
	{
	case T_LINEAR:
		break ;
	case T_LOG:
		d_value = log(d_value) ;
		break ;
	case T_LOG_PLUS1:
		d_value = log(d_value+1) ;
		break ;
	case T_LOG_PLUS3:
		d_value = log(d_value+3) ;
		break ;
	case T_LINEAR_PLUS3:
		d_value = d_value+3 ;
		break ;
	default:
		CIO_DP::message(M_ERROR_DP, "unknown transform\n") ;
		break ;
	}

	INT idx = 0 ;
	REAL ret ;
	for (INT i=0; i<PEN->len; i++)
		if (PEN->limits[i]<=d_value)
			idx++ ;
	
	if (idx==0)
		ret=PEN->penalties[0] ;
	else if (idx==PEN->len)
		ret=PEN->penalties[PEN->len-1] ;
	else
	{
		ret = (PEN->penalties[idx]*(d_value-PEN->limits[idx-1]) + PEN->penalties[idx-1]*
			   (PEN->limits[idx]-d_value)) / (PEN->limits[idx]-PEN->limits[idx-1]) ;  
	}
	//if (p_value>=30 && p_value<150)
	//fprintf(stderr, "%s %i(%i) -> %1.2f\n", PEN->name, p_value, idx, ret) ;
	
	if (PEN->next_pen && follow_next)
		ret+=lookup_penalty(PEN->next_pen, p_value, svm_values);

	return ret ;
}
