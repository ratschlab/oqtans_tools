#ifndef _FILLMATRIX_H__
#define _FILLMATRIX_H__

#include "Mathmatics_dp.h"
#include "penalty_info_dp.h"

#include <algorithm> //Eigentlich in Header
#include <assert.h>
#include <ctype.h>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <string.h>

#include <vector>

struct prev_score { //24B
  double value; //8
  int prev_i; //4
  int prev_j; //4
  int prev_matrix_no;//4
  int num_matches;
  int num_mismatches;
  int num_gaps;
};
typedef struct prev_score Prev_score;

struct position_score{
  int read_pos;
  int dna_pos;
  struct seed_element *next_seed;
  int path_number; //In case of next_seed not NULL
  int path_number_matrices; 
  double partial_score;
  int num_mm;
  int num_gaps;
  int num_introns;
};
typedef struct position_score PosScore;

//Structure for an alignment from a seed
struct seed_element{
	int read_pos;
	int dna_pos;
	int max_gaps;
	int max_mm;
	int max_introns;
	struct prev_score** matrices;
	double* best_scores;
	double best_prev_score;
	PosScore** best_score_pos;
};
typedef struct seed_element SeedElem;


struct splice_pos{
  
  int site;
  int matrix_pos;
  int i;
  int number;
  int* matrices;
};


//extern void print_align(Pre_score* matrix, int length_est,  int length_dna, Align_pair* vektor, int result_length, int print_matrix);

void clean_seed_matrix_vector(std::vector<SeedElem *>& matrix, int nr_paths);
void fast_fill_matrix(int nr_paths_par, int*max_score_positions, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, 
		      double* matchmatrix, penalty_struct* qualityScores, double* donor, double* acceptor, bool remove_duplicate_scores,int seed_i, int seed_j, 
		       std::vector<SeedElem *>& seed_matrix_left, std::vector<SeedElem *>& seed_matrix_right, int max_number_introns, 
					  int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity,mode currentMode, bool remapping);


int check_char(char base);

double getScore(double *matchmatrix,int mlen, int dnaChar, int estChar);
inline double getScore(struct penalty_struct* qualityScores, int mlen, int dnaChar, int estChar, double baseScore) 
{
	//printf("Starting scoring / scheme is USE_QUALITY_SCORES\n");
	//printf("estChar,dnaChar are: %d,%d",estChar,dnaChar);

	//assert(estChar > 0 && estChar < 6);
	//assert(dnaChar > -1 && dnaChar < 6);
   
   int currentPos = (estChar-1)*6+dnaChar;
   struct penalty_struct *currentPen = &qualityScores[currentPos];
   //fprintf(stderr, "currentPos=%i %ld\n", currentPos, currentPen) ;
   
   double score = 0 ;

   if ( currentPen->cache!=NULL )
   {
	   score = currentPen->cache[(int)baseScore] ;
	   //assert( fabs(baseScore-floor(baseScore))<1e-5 && baseScore>=0 ) ;
	   return score ;
   }

   //if( estChar == 1 and dnaChar == 1)
   //   printf("A-A : %f %f\n",currentPen.penalties[0],currentPen.penalties[1]);

	INT idx = 0 ;
	for (INT i=0; i<currentPen->len; i++){
	  //fprintf(stdout,"currentPen->len %i last element %f\n",currentPen->len,currentPen->limits[currentPen->len-1]);
	  if (currentPen->limits[i] <= baseScore)
	    idx++ ;
	}
	if (idx==0)
		score=currentPen->penalties[0] ;
	else if (idx==currentPen->len)
		score=currentPen->penalties[currentPen->len-1] ;
	else
	{
		score = (currentPen->penalties[idx]*(baseScore-currentPen->limits[idx-1]) + currentPen->penalties[idx-1]*
			   (currentPen->limits[idx]-baseScore)) / (currentPen->limits[idx]-currentPen->limits[idx-1]) ;  
   }
   return score;
}


extern int char_map[133] ;
inline int check_char(char base)
{
	return char_map[(int)base] ;
}

#endif // _FILLMATRIX_H__
