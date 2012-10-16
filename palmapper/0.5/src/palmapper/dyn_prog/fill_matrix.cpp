// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

/*********************************************************************************************/
// fast_fill_matrix 
// Not fills the entire matrix: based on a banded semi-global alignment algorithm adapted to spliced alignment
// Applies twice the algorithm from a "seed position" computed from a long hit of GenomeMapper
// First unspliced alignment and in unspliced alignment, first in the diagonal from seed position 
/*********************************************************************************************/

/*
  matchmatrix: columnwise
  - A C G T N (dna)
  - x x x x x x
  A x x x x x x	
  C x x x x x x
  G x x x x x x
  T x x x x x x
  N x x x x x x
  (est)
*/

/*
  alignment matrix: columnwise
  |DNA . . .
  -+------------
  R|00 ...
  E|0.
  A|. .
  D|.  .
  .|.
  .|
  .|
*/


int number_fill_matrix;

static const int MAX_SPLICE_MISMATCH_QUAL_SINGLE=30 ;
static const int MAX_SPLICE_MISMATCH_QUAL_TOTAL=60 ;
static const int MAX_SPLICE_MISMATCH_QUAL_EXTEND=2 ;
static const int MAX_SPLICE_MISMATCH_NUM_N=2 ;

#include "fill_matrix.h"
#include "debug_tools.h"

#define D_USE_QUALITY_SCORES 1

inline bool myisfinite(double x)
{
	if (x<=-ALMOST_INFINITY || x>=ALMOST_INFINITY)
		return false ;
	return true ;
}

inline bool isnotminusinf(double x)
{
	if (x<=-ALMOST_INFINITY)
		return false ;
	return true ;
}

int char_map[133]={-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, 5, -1, -1, -1, -1, -1, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1} ;


/* system of Palma. Otherwise scoring_functions should point to an array of
 * plifs scoring the respective pairs of characters together with the EST
 * quality score.
 */





void initScoreCache(struct penalty_struct* qualityScores, int mlen) 
{
	const int max_qual = 100 ;
	for (int estChar=1; estChar<6; estChar++)
		for (int dnaChar=0; dnaChar<6; dnaChar++)
		{
			int currentPos = (estChar-1)*6+dnaChar;
			struct penalty_struct *currentPen = &qualityScores[currentPos];
			if (currentPen->cache)
				continue ;

			double * cache = (double*)malloc(sizeof(double)*max_qual) ;
			for (int qual=0; qual<max_qual; qual++)
				cache[qual]=getScore(qualityScores, mlen, dnaChar, estChar, qual) ;

			currentPen->cache = cache ;
		}
}


void clean_seed_matrix_vector(std::vector<SeedElem*> &matrix, int nr_paths){

	//  fprintf(stdout,"Clean seed matrix...\n");
	for (int n = 0; n < (int)matrix.size(); n++){
		if (matrix[n]!=NULL){
			for (int p=nr_paths-1;p>=0;p--){
				delete[] matrix[n]->matrices[p];
				if (matrix[n]->best_score_pos[p]!=NULL)
					matrix[n]->best_score_pos[p]->next_seed=NULL;
				delete matrix[n]->best_score_pos[p];
			}
			delete[] matrix[n]->matrices;
			delete[] matrix[n]->best_score_pos;
			delete[] matrix[n]->best_scores;
			delete matrix[n];
		} 
	}
	matrix.clear();
	// fprintf(stdout,"Clean seed matrix...END\n");
}


// Sort best_scores and best_score_pos for the input seed_matrix
void sort_best_scores(SeedElem* seed, int nr_paths){

	double last_score=seed->best_scores[nr_paths-1];
  
	PosScore* temp_posscore=NULL;
	int z=nr_paths-2;

	while(z>=0 && last_score>seed->best_scores[z]){
    
		temp_posscore=seed->best_score_pos[z+1];
		seed->best_score_pos[z+1]=seed->best_score_pos[z];
		seed->best_score_pos[z]=temp_posscore;
		seed->best_scores[z+1]=seed->best_scores[z];
		seed->best_scores[z]=last_score;
		z--;
	}
	temp_posscore=NULL;

}


// Sort paths for a given position in the matrix according to scores
void sort_best_paths(Prev_score*matrices[], int nr_paths,int matrix_position){

	double last_score=((Prev_score*)matrices[nr_paths-1]+matrix_position)->value;
  
	Prev_score temp_prevscore;
	int z=nr_paths-2;

	while(z>=0 && last_score>((Prev_score*)matrices[z]+matrix_position)->value){
		temp_prevscore.value=last_score;
		temp_prevscore.prev_i=((Prev_score*)matrices[z+1]+matrix_position)->prev_i;
		temp_prevscore.prev_j=((Prev_score*)matrices[z+1]+matrix_position)->prev_j;
		temp_prevscore.prev_matrix_no=((Prev_score*)matrices[z+1]+matrix_position)->prev_matrix_no;
		temp_prevscore.num_matches=((Prev_score*)matrices[z+1]+matrix_position)->num_matches;
		temp_prevscore.num_mismatches=((Prev_score*)matrices[z+1]+matrix_position)->num_mismatches;
		temp_prevscore.num_gaps=((Prev_score*)matrices[z+1]+matrix_position)->num_gaps;
    
		((Prev_score*)matrices[z+1]+matrix_position)->value= ((Prev_score*)matrices[z]+matrix_position)->value;
		((Prev_score*)matrices[z+1]+matrix_position)->prev_i= ((Prev_score*)matrices[z]+matrix_position)->prev_i;
		((Prev_score*)matrices[z+1]+matrix_position)->prev_j= ((Prev_score*)matrices[z]+matrix_position)->prev_j;
		((Prev_score*)matrices[z+1]+matrix_position)->prev_matrix_no=((Prev_score*)matrices[z]+matrix_position)->prev_matrix_no;
		((Prev_score*)matrices[z+1]+matrix_position)->num_matches=((Prev_score*)matrices[z]+matrix_position)->num_matches;
		((Prev_score*)matrices[z+1]+matrix_position)->num_mismatches=((Prev_score*)matrices[z]+matrix_position)->num_mismatches;
		((Prev_score*)matrices[z+1]+matrix_position)->num_gaps=((Prev_score*)matrices[z]+matrix_position)->num_gaps;
    
    
		((Prev_score*)matrices[z]+matrix_position)->value=temp_prevscore.value;
		((Prev_score*)matrices[z]+matrix_position)->prev_i=temp_prevscore.prev_i;
		((Prev_score*)matrices[z]+matrix_position)->prev_j=temp_prevscore.prev_j;
		((Prev_score*)matrices[z]+matrix_position)->prev_matrix_no=temp_prevscore.prev_matrix_no;
		((Prev_score*)matrices[z]+matrix_position)->num_matches=temp_prevscore.num_matches;
		((Prev_score*)matrices[z]+matrix_position)->num_mismatches=temp_prevscore.num_mismatches;
		((Prev_score*)matrices[z]+matrix_position)->num_gaps=temp_prevscore.num_gaps;
    
		z--;
	}
}


int check_min_matches(SeedElem* seed, int nr_paths, int matrix_position, int min_matches, int*matrices, int i, int j, int prev_shift, char* read, int read_len, char* dna, int dna_len, double* read_scores, int verbosity)
{
	double prevValue;
	int num;
	int matrix_number=0;
  

  
	for (int z=0;z<nr_paths;z++)
    {
		int diff_i = min_matches ;

		prevValue = ((Prev_score*)seed->matrices[z] +matrix_position)->value ; 	    
		num = ((Prev_score*)seed->matrices[z] +matrix_position)->num_matches;

		bool conserved_seq=true ;
		int num_N=0 ;
		int num_mismatch_score=0 ;
		int num_mismatch_extend=0 ;
		int pos_i=i+prev_shift;
		int pos_j=j+prev_shift;
		for (int mn=1; mn<=diff_i; mn++)
		{
			// fprintf(stdout,"[check_min_matches] From position %i-%i with %c-%c\n",pos_i,pos_j,read[pos_i],dna[pos_j]);
			if (pos_i<0 || pos_j<0 || pos_i>=read_len || pos_j>=dna_len)
			{
				//fprintf(stdout,"[check_min_matches] out of bounds\n");
				conserved_seq=false ;
				break ;
			}
			if (read[pos_i]!=dna[pos_j])
			{
				if (check_char(read[pos_i])!=5 && check_char(dna[pos_j])!=5 && read_scores[pos_i]>=MAX_SPLICE_MISMATCH_QUAL_SINGLE)
				{
					conserved_seq=false ;
					//fprintf(stdout,"[check_min_matches] mismatch with high quality\n");
					break ;
				}
				else
				{
					diff_i++ ;
					if (check_char(read[pos_i])==5 || check_char(dna[pos_j])==5)
						num_N++ ;
					else
					{
						num_mismatch_score+=read_scores[pos_i] ;
						num_mismatch_extend++ ;
					}
					if (num_N>MAX_SPLICE_MISMATCH_NUM_N || num_mismatch_score>MAX_SPLICE_MISMATCH_QUAL_TOTAL || num_mismatch_extend>MAX_SPLICE_MISMATCH_QUAL_EXTEND) 
					{
						//fprintf(stdout,"[check_min_matches] too many N or too many mismatches or num_mismatch_score too high\n");
						conserved_seq=false;
						break;  
					}
				}
			}
			pos_i+=prev_shift ;
			pos_j+=prev_shift ;
		}

		if (!conserved_seq && num>=min_matches+num_N && num_N<=MAX_SPLICE_MISMATCH_NUM_N){
			if (verbosity>0)
				fprintf(stdout,"[check_min_matches] PROBLEM with number of matches from position %i-%i (%i for %i) with %i N\n",i,j,num,min_matches,num_N);
		}
		if (isnotminusinf(prevValue) && conserved_seq)
		{
			//fprintf(stdout,"Possible matrix: %i with value %f\n",z,prevValue);
			matrices[matrix_number]=z;
			matrix_number++;
		}
    }
	return matrix_number;
}

void print_restricted_matrix(Prev_score* matrices[],int nr_paths, int matrix_len, int row_len){

	for(int z=0; z<nr_paths;z++){
		fprintf(stdout,"Matrix %i:\n",z);
		for (int i=0; i<matrix_len;i++){
			fprintf(stdout,"%f(%i,%i,%i) ", ((Prev_score*)matrices[z]+i)->value,((Prev_score*)matrices[z]+i)->num_gaps,((Prev_score*)matrices[z]+i)->num_mismatches,((Prev_score*)matrices[z]+i)->num_matches);
			if ((i+1)%row_len==0)
				fprintf(stdout,"\n");
      
		}
	}
	fprintf(stdout,"\n\n");
}



void fast_fill_side_unspliced_first(int nr_paths_par,  std::vector<SeedElem*> &seed_matrix, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, double* matchmatrix, penalty_struct* qualityScores, double* main_site_scores, double* comp_site_scores, std::vector<int>& comp_sites,	int seed_read, int seed_dna, double* best_match_scores, bool right_side,bool first_seed,int max_number_introns,	int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity, mode currentMode, bool remapping)
{

//fprintf(stdout,"START: Fill %s side of the matrix from position %i-%i (%i,%i,%i,num_intron=%i)...\n",(right_side==true)?"right":"left",seed_read, seed_dna,max_gap,max_mism,max_edit_op, max_number_introns);
  
	/***************************************************/
	/*Initialization */
	/***************************************************/
	const int mlen = 6; // length of matchmatrix
	int readChar, dnaChar; // 0:'-', 1:'a', 2:'c', 3:'g', 4:'t', 5:'n'
  
	double baseScore=0;
	double *read_scores = prb ;

	double prevValue;
	double tempValue;
	double putativeValue;
	double globalValue;
	int prevGaps;
	int prevMism;

	int matrix_position;
	int matrix_prev_position;
  
	int prev_shift;
	if(right_side)
		prev_shift=-1;
	else
		prev_shift=1;

	// Boolean that indicates if complementary splice sites had been scanned before
	bool comp_sites_filled=false;
  
	// Possible paths with enough matches before a given splice site
	int * possible_matrices= new int[nr_paths_par];
	// Possible splice sites for starting a new alignment
	std::vector<splice_pos *> possible_sites;
  
	if (currentMode==USE_QUALITY_SCORES)
		assert((int)functions->limits[functions->len-1]==2);
  
	//width of the band according to the number of gaps
	int row_len=max_gap*2+1;

	//Length of read to align from seed position and of dna from seed_position to beginning or end
	int i_len;
	int j_len;
	if (right_side){
		i_len=read_len-seed_read;
		j_len=dna_len-seed_dna;
	}
	else{
		i_len=seed_read+1;
		j_len=seed_dna+1;
	}
	int matrix_len=row_len*i_len;

	//Init score matrices
	Prev_score** matrices = new Prev_score*[nr_paths_par];
	for (int z=0; z<nr_paths_par; z++) {
		matrices[z] = new Prev_score[matrix_len]; //Allocation only for the band around the seed position
		memset(matrices[z],0,sizeof(Prev_score)*matrix_len);
	}

	// Diagonal activation: first, only activating of diagonal from seed position
	bool* disabled_diagonal=new bool[row_len];
	for (int d=0;d<max_gap;d++){
		disabled_diagonal[d]=true;
		disabled_diagonal[row_len-1-d]=true;
	}
	disabled_diagonal[max_gap]=false;

	//Init seed position score in matrices
	//Best matrix
	dnaChar=check_char(dna[seed_dna]) ;
	readChar = check_char(read[seed_read]) ;
  
	if (currentMode == USE_QUALITY_SCORES){
		baseScore = read_scores[seed_read];
		((Prev_score*)matrices[0]+ max_gap)->value = getScore(qualityScores,mlen,dnaChar,readChar,baseScore);
	}
	else{
		((Prev_score*)matrices[0]+ max_gap)->value = (matchmatrix[mlen* dnaChar +readChar]);
	}
	
	((Prev_score*)matrices[0]+ max_gap)->prev_i = seed_read+prev_shift;
	((Prev_score*)matrices[0]+ max_gap)->prev_j = seed_dna+prev_shift;
	((Prev_score*)matrices[0]+ max_gap)->prev_matrix_no = 0;
	((Prev_score*)matrices[0]+ max_gap)->num_matches = ((read[seed_read]==dna[seed_dna] || check_char(read[seed_read])==5 || check_char(dna[seed_dna])==5) ? 1:0);
	((Prev_score*)matrices[0]+ max_gap)->num_mismatches = 1-((Prev_score*)matrices[0]+ max_gap)->num_matches ;
	((Prev_score*)matrices[0]+ max_gap)->num_gaps =0;
  
	//Other matrices: score at -INF
	for(int z=1; z<nr_paths_par;z++){
		((Prev_score*)matrices[z]+ max_gap)->value = log(0.0); // -inf
		((Prev_score*)matrices[z]+ max_gap)->prev_i = 0; //seed_read+prev_shift;
		((Prev_score*)matrices[z]+ max_gap)->prev_j = 0; //seed_dna+prev_shift;
		((Prev_score*)matrices[z]+ max_gap)->prev_matrix_no = 0;
		((Prev_score*)matrices[z]+ max_gap)->num_matches = 0;
		((Prev_score*)matrices[z]+ max_gap)->num_mismatches = 0;
		((Prev_score*)matrices[z]+ max_gap)->num_gaps = 0;
	}

	//Alignment from the seed with (seed_read, seed_dna) position
	SeedElem* current_seed= new SeedElem();
	current_seed->read_pos=seed_read;
	current_seed->dna_pos=seed_dna;
	current_seed->max_gaps=max_gap;
	current_seed->max_mm=max_mism;
	current_seed->max_introns=max_number_introns;
	current_seed->best_scores= new double[nr_paths_par];
	current_seed->best_prev_score= log(0.0);
	current_seed->best_score_pos=new PosScore*[nr_paths_par];
	for(int z =0; z<nr_paths_par;z++){
		current_seed->best_scores[z]=log(0.0);
		PosScore * pscore= new PosScore();    
		pscore->read_pos=0;
		pscore->dna_pos=0;
		pscore->num_gaps=0;
		pscore->num_mm=0;
		pscore->num_introns=0;
		pscore->next_seed=NULL;
		pscore->path_number=0; //path number in next matrix
		pscore->path_number_matrices=0;
		pscore->partial_score=0; //score from seed position to (read_pos,dna_pos)
		current_seed->best_score_pos[z]=pscore;
	}
	current_seed->matrices=matrices;
	seed_matrix.push_back(current_seed);
  

	//Can align all read nucleotides: all read nucleotides are covered by the band
	if (j_len+max_gap >= i_len){
    
		/***************************************************/
		/*banded SW algorithm                              */
		/***************************************************/
		for(int ni=0; ni<i_len;ni++){
      
			//Position in the read
			int i= seed_read-prev_shift*ni;

			//dna interval for this row
			int left_bound =seed_dna+(i-seed_read)-max_gap;
			int right_bound=seed_dna+(i-seed_read)+max_gap;
      
			for (int nj=0; nj<row_len;nj++){
				int j;
				if (right_side)
					j=left_bound+nj;
				else
					j=right_bound-nj;

				matrix_position= ni*row_len+nj;


				//Diagonal not active or position out of bounds
				if(disabled_diagonal[nj] || (right_side && (j<seed_dna || j>=dna_len))||(!right_side &&(j>seed_dna || j<0))){
					for(int z=0;z<nr_paths_par;z++){
						((Prev_score*)matrices[z] + matrix_position)->value = log(0.0) ;
						((Prev_score*)matrices[z] + matrix_position)->prev_i = 0;
						((Prev_score*)matrices[z] + matrix_position)->prev_j = 0;
						((Prev_score*)matrices[z] + matrix_position)->prev_matrix_no = 0;
						((Prev_score*)matrices[z] + matrix_position)->num_matches = 0;
						((Prev_score*)matrices[z] + matrix_position)->num_mismatches = 0;
						((Prev_score*)matrices[z] + matrix_position)->num_gaps = 0;
					}
					continue;
				}


				// seed position filled before
				if(i==seed_read && j==seed_dna){
					if(ni==i_len-1){
						current_seed->best_scores[0]=((Prev_score*)current_seed->matrices[0] +matrix_position)->value; 	    
						current_seed->best_score_pos[0]->read_pos=i;
						current_seed->best_score_pos[0]->dna_pos=j;
						current_seed->best_score_pos[0]->num_gaps=0;
						current_seed->best_score_pos[0]->num_mm=((Prev_score*)current_seed->matrices[0] +matrix_position)->num_mismatches;
						current_seed->best_score_pos[0]->num_introns=0;
						current_seed->best_score_pos[0]->next_seed=NULL;
						current_seed->best_score_pos[0]->path_number=0;
						current_seed->best_score_pos[0]->path_number_matrices=0;
						current_seed->best_score_pos[0]->partial_score=((Prev_score*)current_seed->matrices[0] +matrix_position)->value;
					}
					continue;
				}

	
				/********************************************************/
				/* 1. DETECTION OF POSSIBLE SPLICE SITE                 */
				/********************************************************/
				if (main_site_scores[j]>-ALMOST_INFINITY && main_site_scores[j]<ALMOST_INFINITY){	
					if ((max_number_introns>0) && (matrix_position-row_len>=0)){

						//Compute minimum number of matches needed with relaxing the constraint close to the first seed position
						int diff_i;
						if (first_seed){
							diff_i=prev_shift*(seed_read-i);
							if (diff_i>min_match)
								diff_i=min_match;
						}
						else
							diff_i=min_match;

						//Number of paths that has at least diff_i consecutive matches that lead to the previous position in the diagonal
						int number;
						number=check_min_matches(current_seed,nr_paths_par,matrix_position-row_len, diff_i,possible_matrices, i, j, prev_shift, read, read_len, dna, dna_len, read_scores, verbosity);
						if (number>0){
							splice_pos* sp= new splice_pos();
							sp->site=j; //position of the splice site
							sp->i=i+prev_shift; //position of the previous position in the diagonal
							sp->matrix_pos=matrix_position-row_len; //position in the matrix of the previous position in the diagonal
							sp->number=number;
							sp->matrices=new int[number];
							for(int e=0;e<number;e++)
								sp->matrices[e]=possible_matrices[e];
							possible_sites.push_back(sp);
						}
					}
				}
				/*1. END */
	


				/*************************************************************/
				/* 2. FILLING THE DIAGONAL (MATCH/MISMATCH)                  */
				/* right side: (i-1,j-1)->(i,j); left side: (i+1,j+1)->(i,j) */
				/*************************************************************/
				dnaChar = check_char(dna[j]);
				readChar = check_char(read[i]);
				assert(dnaChar!=-1 && readChar!=-1);

				if (currentMode == USE_QUALITY_SCORES)
					baseScore = read_scores[i];
				
				// Best score of what it leaves to align
				if (ni<i_len-1){
					if(right_side)
						putativeValue= best_match_scores[i+1];
					else
						putativeValue= best_match_scores[0]-best_match_scores[i];
				}
				else
					putativeValue=0;

				//Least good best global score known
				globalValue= current_seed->best_scores[nr_paths_par-1];
	
				//Information about previous position in the diagonal
				matrix_prev_position= matrix_position-row_len;
				int num_unfilled=0;
				for(int z=0;z<nr_paths_par;z++){
					Prev_score* actMatrix = (Prev_score*)matrices[z]; 
					prevValue = ((Prev_score*)actMatrix +matrix_prev_position)->value ; 
					prevGaps=((Prev_score*)actMatrix +matrix_prev_position)->num_gaps;
					prevMism=((Prev_score*)actMatrix +matrix_prev_position)->num_mismatches;
	
					if (isnotminusinf(prevValue)){

						if (currentMode == USE_QUALITY_SCORES)
							tempValue = prevValue + getScore(qualityScores,mlen,dnaChar,readChar,baseScore);
						else
							tempValue = prevValue +(matchmatrix[mlen* dnaChar +readChar]);

						//Fill if tempValue is greater of equal to an existing spliced alignment
						//If mismatch, does not have to rise above the number of allowed mismatches and edit operations 
						if (isnotminusinf(tempValue) && (read[i]==dna[j] || (prevMism+1<=max_mism && prevMism+prevGaps+1<=max_edit_op)) && tempValue+putativeValue>globalValue){
							((Prev_score*)actMatrix + matrix_position)->value = tempValue;
							((Prev_score*)actMatrix + matrix_position)->prev_i = i+prev_shift; 
							((Prev_score*)actMatrix + matrix_position)->prev_j = j+prev_shift; 
							((Prev_score*)actMatrix + matrix_position)->prev_matrix_no = z;
							if (read[i]==dna[j] || check_char(read[i])==5 || check_char(dna[j])==5) {
								((Prev_score*)actMatrix + matrix_position)->num_matches = ((Prev_score*)actMatrix +matrix_prev_position)->num_matches+1; 
								((Prev_score*)actMatrix + matrix_position)->num_mismatches = prevMism;
							}
							else{
								((Prev_score*)actMatrix + matrix_position)->num_matches = 0 ;
								((Prev_score*)actMatrix + matrix_position)->num_mismatches = prevMism+1;
							}	      
							((Prev_score*)actMatrix + matrix_position)->num_gaps = prevGaps;

							//Last row to fill: check now if unspliced alignment is better
							if(ni==i_len-1){
								if (tempValue>current_seed->best_scores[nr_paths_par-1]){
									current_seed->best_scores[nr_paths_par-1]=tempValue;
									current_seed->best_score_pos[nr_paths_par-1]->read_pos=i;
									current_seed->best_score_pos[nr_paths_par-1]->dna_pos=j;
									current_seed->best_score_pos[nr_paths_par-1]->num_gaps=((Prev_score*)actMatrix + matrix_position)->num_gaps;									
									current_seed->best_score_pos[nr_paths_par-1]->num_mm=((Prev_score*)actMatrix + matrix_position)->num_mismatches;
									current_seed->best_score_pos[nr_paths_par-1]->num_introns=0;
									current_seed->best_score_pos[nr_paths_par-1]->next_seed=NULL;
									current_seed->best_score_pos[nr_paths_par-1]->path_number=0;
									current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=z;
									current_seed->best_score_pos[nr_paths_par-1]->partial_score=tempValue;
									// Resort best_score_pos and best_scores
									sort_best_scores(current_seed, nr_paths_par);
								}
							}
						}
						//Too much mismatches/edit operations or better alignment with intron: desactivation of the diagonal
						else{
							((Prev_score*)actMatrix + matrix_position)->value = log(0.0) ;
							((Prev_score*)actMatrix + matrix_position)->prev_i = 0;
							((Prev_score*)actMatrix + matrix_position)->prev_j = 0;
							((Prev_score*)actMatrix + matrix_position)->prev_matrix_no = 0;
							((Prev_score*)actMatrix + matrix_position)->num_matches = 0;
							((Prev_score*)actMatrix + matrix_position)->num_mismatches = 0;
							((Prev_score*)actMatrix + matrix_position)->num_gaps = 0;

							num_unfilled++;
						}
					}
					else{
						((Prev_score*)actMatrix + matrix_position)->value = log(0.0) ;
						((Prev_score*)actMatrix + matrix_position)->prev_i = 0;
						((Prev_score*)actMatrix + matrix_position)->prev_j = 0;
						((Prev_score*)actMatrix + matrix_position)->prev_matrix_no = 0;
						((Prev_score*)actMatrix + matrix_position)->num_matches = 0;
						((Prev_score*)actMatrix + matrix_position)->num_mismatches = 0;
						((Prev_score*)actMatrix + matrix_position)->num_gaps = 0;
	    
						num_unfilled++;
					}
				}
				// If no solution for all paths, desactivate this diagonal
				if(num_unfilled==nr_paths_par)
					disabled_diagonal[nj]=true;
				/*2. END */


				/******************************************************************************************/
				/* 3. POSITION FILLED WITH A MISMATCH OR NOT FILLED                                       */
				/* previous position in the diagonal can lead to a gap on read or dna instead of mismatch */
				/* Also look at possible splice sites stored from this diagonal                           */
				/******************************************************************************************/
				if (disabled_diagonal[nj] || read[i]!=dna[j]){
					//if (possible_sites.size()>10)
					//fprintf(stdout, "possible_sites.size()=%i\n", (int)possible_sites.size()) ;
					/* A. Explore splice sites */

					if(!remapping){
						for(int ss=(int)possible_sites.size()-1;ss>=0;ss--){
	    
							//Position of the splice site 'G' of 'GT/C' or 'G' of 'AG' on DNA
							int posj=((splice_pos*)possible_sites[ss])->site;

							//Last position before splice site
							int posi=possible_sites[ss]->i;
							int matrix_pos=((splice_pos*)possible_sites[ss])->matrix_pos;

							int number=possible_sites[ss]->number;
							int* pmatrices=possible_sites[ss]->matrices;
	    
							
							if (matrix_pos%row_len==nj){

								prevGaps=((Prev_score*)matrices[0] + matrix_pos)->num_gaps;
								prevMism=((Prev_score*)matrices[0] + matrix_pos)->num_mismatches;
	      
	      
								//Figure out possible complementary splice sites if not done before
								if (first_seed && !comp_sites_filled){		
									comp_sites_filled=true;
									int first_ss=((splice_pos*)possible_sites[0])->site;
									if(right_side){
										for(int comp_ss=first_ss+2;comp_ss<dna_len-1;comp_ss++){
											if (comp_site_scores[comp_ss]>-ALMOST_INFINITY && comp_site_scores[comp_ss]<ALMOST_INFINITY)
												comp_sites.push_back(comp_ss);
										}
									}
									else{
										for(int comp_ss=first_ss-2;comp_ss>0;comp_ss--){
											if (comp_site_scores[comp_ss]>-ALMOST_INFINITY && comp_site_scores[comp_ss]<ALMOST_INFINITY)
												comp_sites.push_back(comp_ss);
										}
									}
								}

	      
								//Alignment from each complementary splice site
								for(int comp_ss=0; comp_ss<(int)comp_sites.size();comp_ss++){
	      
									// Next complementary splice site relative to the current splice site
									int jj=comp_sites[comp_ss]; 
									if (-prev_shift*jj<=posj*-prev_shift)
										continue;

									//Compute minimum number of matches needed with relaxing the constraint at the ends of the read
									int diff_i;
									if(right_side)
									{
										diff_i=read_len-1-posi;
										if (posi==read_len-1 || diff_i>min_match)
											diff_i=min_match;
									}
									else
									{
										diff_i=posi;
										if (posi==0 || posi>min_match)
											diff_i=min_match;
									}
		
									//Check at least diff_i consecutive matches after this complementary splice site
									bool conserved_seq=true;	      
									int conserved_seq_mismatches=0 ;
									int ii=posi;
									int num_N=0 ;
									int num_mismatch_score=0 ;
									int num_mismatch_extend=0 ;
									for(int nm=1; nm<=diff_i;nm++)
									{
										ii=ii-prev_shift; 
										jj=jj-prev_shift; // One position after 'G' of 'AG' or before 'G' of 'GT/C'
										if ( (right_side && (ii>=read_len || jj>=dna_len)) || (!right_side && (ii<0 || jj<0)))
										{
											conserved_seq_mismatches++ ;
											conserved_seq=false;
											break;		   
										}
										else{
											if (read[ii]!=dna[jj]) 
											{
												if (check_char(read[ii])!=5 && check_char(dna[jj])!=5 && read_scores[ii]>=MAX_SPLICE_MISMATCH_QUAL_SINGLE)
												{
													conserved_seq=false;
													break;
												}
												else
												{
													diff_i++ ;
													if (check_char(read[ii])==5 || check_char(dna[jj])==5)
														num_N++ ;
													else
													{
														num_mismatch_score+=read_scores[ii] ;
														num_mismatch_extend++ ;
													}
													if (num_N>MAX_SPLICE_MISMATCH_NUM_N || num_mismatch_score>MAX_SPLICE_MISMATCH_QUAL_TOTAL || num_mismatch_extend>MAX_SPLICE_MISMATCH_QUAL_EXTEND)
													{
														conserved_seq=false;
														break;  
													}
												}
											}
										}
									}

									if (conserved_seq){
		  


										int seed_already_filled=-1;
										bool continue_searching=true;
										double tempSplicedScore = main_site_scores[posj]+comp_site_scores[comp_sites[comp_ss]];
										double currentScore= ((Prev_score*)current_seed->matrices[pmatrices[0]] +matrix_pos)->value ; 	    

										for(int num_seed=0; num_seed< (int)seed_matrix.size();num_seed++){
										
											if(seed_matrix[num_seed]->read_pos==posi-prev_shift && 
											   seed_matrix[num_seed]->dna_pos==comp_sites[comp_ss]-prev_shift && seed_matrix[num_seed]->max_introns==max_number_introns-1)
											{
								
												if ((seed_matrix[num_seed]->best_score_pos[0]->num_gaps>max_gap-prevGaps && seed_matrix[num_seed]->best_score_pos[0]->num_mm>=max_mism-prevMism)||
													(seed_matrix[num_seed]->best_score_pos[0]->num_gaps>=max_gap-prevGaps && seed_matrix[num_seed]->best_score_pos[0]->num_mm>max_mism-prevMism)){
													continue_searching=false;
													continue;												
												}
											
												if ((seed_matrix[num_seed]->max_gaps>=max_gap-prevGaps && seed_matrix[num_seed]->max_mm>=max_mism-prevMism)&&
													(seed_matrix[num_seed]->best_score_pos[0]->num_gaps<=max_gap-prevGaps && seed_matrix[num_seed]->best_score_pos[0]->num_mm<=max_mism-prevMism)){
												
													if (seed_matrix[num_seed]->best_scores[0] > log(0.0) && tempSplicedScore+currentScore >= seed_matrix[num_seed]->best_prev_score){
														continue_searching=true;
														seed_already_filled=num_seed;
														break;
													}
													else{													
														continue_searching=false;
														break;
													}												
												}

												if (seed_matrix[num_seed]->max_gaps < max_gap-prevGaps && seed_matrix[num_seed]->max_mm < max_mism-prevMism){		
													if (tempSplicedScore+currentScore >= seed_matrix[num_seed]->best_prev_score){
														continue_searching=true;
														continue;
													}
													else{													
														continue_searching=false;
														break;
													}												
												}
											}
										}
									
														  
										// if (!continue_searching)									   
										// 	fprintf(stdout,"**1** Don't fill this seed matrix and don't compare results\n");
										// else{
										// 	if (seed_already_filled!=-1)
										// 		fprintf(stdout,"**2** Use already filled matrix\n");
										// 	//else
										// 		//fprintf(stdout,"**3** New seed matrix\n");
										// }
									
										if (continue_searching){
										
											if(seed_already_filled==-1){
												//Number of this seedElem in the vector seedMatrix (because of recursive calls that add new seedElem)
												seed_already_filled=seed_matrix.size();
												fast_fill_side_unspliced_first(nr_paths_par, seed_matrix, read_len, dna_len, read, dna, prb, functions, matchmatrix, qualityScores, 
																			   main_site_scores, comp_site_scores, comp_sites, posi-prev_shift,
																			   comp_sites[comp_ss]-prev_shift, best_match_scores, right_side,false,
																			   max_number_introns-1,max_gap-prevGaps,max_mism-prevMism,max_edit_op-(prevGaps+prevMism),min_match, verbosity,currentMode, remapping);
											}

											//Keep best scores
											//double tempSplicedScore = main_site_scores[posj]+comp_site_scores[comp_sites[comp_ss]];
											//fprintf(stdout,"Intron score %f %f \n",main_site_scores[posj],comp_site_scores[comp_sites[comp_ss]]);
		  
											for (int z=0; z<nr_paths_par;z++){
												for (int zz=0; zz<number;zz++){
													double priorScore= ((Prev_score*)current_seed->matrices[pmatrices[zz]] +matrix_pos)->value ; 	    
													//fprintf(stdout,"prior score %f new score %f splice score %f\n",priorScore,seed_matrix[seed_already_filled]->best_scores[z],tempSplicedScore);
												
													if(tempSplicedScore +priorScore>seed_matrix[seed_already_filled]->best_prev_score)
														seed_matrix[seed_already_filled]->best_prev_score=tempSplicedScore +priorScore;
												
													if(tempSplicedScore +priorScore+seed_matrix[seed_already_filled]->best_scores[z] > current_seed->best_scores[nr_paths_par-1]){											
														current_seed->best_scores[nr_paths_par-1]=tempSplicedScore+priorScore+seed_matrix[seed_already_filled]->best_scores[z];
														current_seed->best_score_pos[nr_paths_par-1]->read_pos=posi;
														current_seed->best_score_pos[nr_paths_par-1]->dna_pos=posj+prev_shift;
														current_seed->best_score_pos[nr_paths_par-1]->num_gaps=((Prev_score*)current_seed->matrices[pmatrices[zz]] +matrix_pos)->num_gaps + seed_matrix[seed_already_filled]->best_score_pos[z]->num_gaps;
														current_seed->best_score_pos[nr_paths_par-1]->num_mm=((Prev_score*)current_seed->matrices[pmatrices[zz]] +matrix_pos)->num_mismatches+ seed_matrix[seed_already_filled]->best_score_pos[z]->num_mm;
														current_seed->best_score_pos[nr_paths_par-1]->num_introns=seed_matrix[seed_already_filled]->best_score_pos[z]->num_introns +1;					
														current_seed->best_score_pos[nr_paths_par-1]->next_seed=seed_matrix[seed_already_filled];
														current_seed->best_score_pos[nr_paths_par-1]->path_number=z;
														current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=pmatrices[zz];
														current_seed->best_score_pos[nr_paths_par-1]->partial_score=priorScore+tempSplicedScore;
														// Resort best_score_pos and best_scores
														sort_best_scores(current_seed, nr_paths_par);
													}
													else{
														break;
													}
												}
											}
										}
									}
								
								}//End look through complementary splice sites

								delete[] possible_sites[ss]->matrices;
								possible_sites[ss]->matrices=NULL;
								delete possible_sites[ss];
								if (ss!=(int)possible_sites.size()-1){
									possible_sites[ss]=possible_sites[possible_sites.size()-1];
								}
								possible_sites.pop_back();
							}
						}
						/* A.END */
					}//remapping
				
	  
					/* B. Allow gaps from the previous position in the diagonal */
					matrix_prev_position= matrix_position-row_len;
	  
					/* B.1. Gap on DNA sequence */
					if (nj!=0 && matrix_prev_position>=0){
						dnaChar = check_char(dna[j+prev_shift]);
						readChar = check_char(read[i]);
						assert(dnaChar!=-1 && readChar!=-1);

						if (currentMode == USE_QUALITY_SCORES)
							baseScore = read_scores[i];
	    
						// Best score of what it leaves to align
						if (ni<i_len-1){
							if(right_side)
								putativeValue= best_match_scores[i+1];
							else
								putativeValue= best_match_scores[0]-best_match_scores[i];
						}
						else
							putativeValue=0;
	    
						//Worst best global score known
						globalValue= current_seed->best_scores[nr_paths_par-1];
						//fprintf(stdout,"Putative value: %f Best score: %f\n",putativeValue,globalValue);
	    

						for(int z=0;z<nr_paths_par;z++){
							Prev_score* actMatrix = (Prev_score*)matrices[z]; 

							//Information about the previous position for path z
							prevGaps=((Prev_score*)actMatrix + matrix_prev_position)->num_gaps;
							prevMism=((Prev_score*)actMatrix + matrix_prev_position)->num_mismatches;
							prevValue = ((Prev_score*)actMatrix +matrix_prev_position)->value ;

							//Gap possible according to the number of gaps and mismatches at matrix_prev_position
							if (prevGaps<max_gap && prevGaps+prevMism<max_edit_op){

								if (currentMode == USE_QUALITY_SCORES)
									tempValue = prevValue + getScore(qualityScores,mlen,0,readChar,baseScore);
								else
									tempValue = prevValue +(matchmatrix[readChar]); /* score(READ,gap) */

								if (isnotminusinf(tempValue)&& tempValue> ((Prev_score*)matrices[nr_paths_par-1] +matrix_position-1)->value && tempValue+putativeValue>globalValue){
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->value = tempValue;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->prev_i = i+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->prev_j = j+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->prev_matrix_no = z;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->num_matches = 0;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->num_mismatches = prevMism;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-1)->num_gaps = prevGaps+1;
									disabled_diagonal[nj-1]=false;
									sort_best_paths(matrices,nr_paths_par,matrix_position-1);

									//Last row to fill: check now if unspliced alignment is better
									if(ni==i_len-1){
										if (tempValue>current_seed->best_scores[nr_paths_par-1]){
											current_seed->best_scores[nr_paths_par-1]=tempValue;
											current_seed->best_score_pos[nr_paths_par-1]->read_pos=i;
											current_seed->best_score_pos[nr_paths_par-1]->dna_pos=j+prev_shift;
											current_seed->best_score_pos[nr_paths_par-1]->num_gaps=((Prev_score*)actMatrix + matrix_position-1)->num_gaps;									
											current_seed->best_score_pos[nr_paths_par-1]->num_mm=((Prev_score*)actMatrix + matrix_position-1)->num_mismatches;
											current_seed->best_score_pos[nr_paths_par-1]->num_introns=0;
											current_seed->best_score_pos[nr_paths_par-1]->next_seed=NULL;
											current_seed->best_score_pos[nr_paths_par-1]->path_number=0;
											current_seed->best_score_pos[nr_paths_par-1]->path_number_matrices=z;
											current_seed->best_score_pos[nr_paths_par-1]->partial_score=tempValue;
											// Resort best_score_pos and best_scores
											sort_best_scores(current_seed, nr_paths_par);
											//worst best global score known
											globalValue= current_seed->best_scores[nr_paths_par-1];
										}
									}
								}
							}
						}
					}/* B.1. END */
	  
					/* B.2. Gap on READ sequence */
					if (nj!=row_len-1 && matrix_prev_position>=0){
						dnaChar = check_char(dna[j]);
						readChar = check_char(read[i+prev_shift]);
						assert(dnaChar!=-1 && readChar!=-1);

						if (currentMode == USE_QUALITY_SCORES)
							baseScore = read_scores[i+prev_shift];

	    
						// Best score of what it leaves to align
						if (ni<i_len){
							if(right_side)
								putativeValue= best_match_scores[i];
							else
								putativeValue= best_match_scores[0]-best_match_scores[i+prev_shift];
						}
						else
							putativeValue=0;

						//Worst best global score known
						globalValue= current_seed->best_scores[nr_paths_par-1];

						for(int z=0;z<nr_paths_par;z++){
							Prev_score* actMatrix = (Prev_score*)matrices[z]; 

							//Information about the previous position for path z
							prevGaps=((Prev_score*)actMatrix + matrix_prev_position)->num_gaps;
							prevMism=((Prev_score*)actMatrix + matrix_prev_position)->num_mismatches;
							prevValue = ((Prev_score*)actMatrix +matrix_prev_position)->value ;

							//Gap possible according to the number of gaps and mismatches at matrix_prev_position
							if (prevGaps<max_gap && prevGaps+prevMism<max_edit_op){
	    
								if (currentMode == USE_QUALITY_SCORES)
									tempValue = prevValue + matchmatrix[dnaChar];
								else
									tempValue = prevValue + matchmatrix[mlen*dnaChar];   /* score(gap,DNA) */

								if (isnotminusinf(tempValue)&& tempValue> ((Prev_score*)matrices[nr_paths_par-1] +matrix_position-row_len+1)->value && tempValue+putativeValue>globalValue){
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->value = tempValue;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->prev_i = i+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->prev_j = j+prev_shift; /* predecessor */
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->prev_matrix_no = z;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->num_matches = 0;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->num_mismatches = prevMism;
									((Prev_score*)matrices[nr_paths_par-1] + matrix_position-row_len+1)->num_gaps = prevGaps+1;
									disabled_diagonal[nj+1]=false;
									sort_best_paths(matrices,nr_paths_par,matrix_position-row_len+1);
								}
							}
						}
					}/* B.2. END */
				}/* 3. END */

			}//end j
	
			bool is_bad_row=true;
			for(int t=0;t<row_len;t++)
				is_bad_row=is_bad_row && (disabled_diagonal[t]==1);
			if (is_bad_row){
				break;
			}
		}//end i
	}
  
	////Print best scores
	//for (int z=0; z<nr_paths_par; z++)
	//   fprintf(stdout,"%i best score: %f\n",z+1,current_seed->best_scores[z]);
	//if(first_seed)
	//   print_restricted_matrix(matrices,nr_paths_par,matrix_len,row_len);
  

	// Free memory
	current_seed=NULL;
	for(int s=0;s<(int)possible_sites.size();s++){
		delete[] possible_sites[s]->matrices;
		possible_sites[s]->matrices=NULL;
		delete possible_sites[s];
	}
	possible_sites.clear();
	delete[] possible_matrices;  
	delete[] disabled_diagonal;

	//fprintf(stdout,"END: Fill %s side of the matrix from position %i-%i (%i,%i,%i,num_intron=%i)...\n",(right_side==true)?"right":"left",seed_read, seed_dna,max_gap,max_mism,max_edit_op, max_number_introns);

	//  fprintf(stdout,"Fill a %s side of the matrix from position %i-%i...END\n",right_side?"right":"left",seed_read, seed_dna);
}



void fast_fill_matrix(int nr_paths_par, int*max_score_positions, int read_len, int dna_len, char* read, char* dna, double* prb, penalty_struct* functions, 
					  double* matchmatrix, penalty_struct* qualityScores, double* donor, double* acceptor, bool remove_duplicate_scores,int seed_i, int seed_j, 
					  std::vector<SeedElem *>& seed_matrix_left, std::vector<SeedElem *>& seed_matrix_right, int max_number_introns, 
					  int max_gap, int max_mism, int max_edit_op, int min_match, int verbosity,mode currentMode, bool remapping)
{
  
	const int MMATRIX_LEN = 6; // length of matchmatrix
  
	//fprintf(stdout,"Max number of exons %i\n",max_number_introns);

	//  printf("Entering fill_matrix...\n");
  
	/*********************************************************************************************/
	/*Best score for a matching alignment from a position i in read sequence   */
	/*********************************************************************************************/
	// fprintf(stdout,"Best match scores...\n");
	double* best_match_scores= new double[read_len];
	double temp_best=0;
	double *read_scores = prb;
  
	for(int i=read_len-1; i >=0;i--){
		
		double score;
		
		if (currentMode == USE_QUALITY_SCORES)
			score=getScore(qualityScores,MMATRIX_LEN,check_char(read[i]),check_char(read[i]),read_scores[i]);
		else
			score=matchmatrix[MMATRIX_LEN*check_char(read[i])+check_char(read[i])];
		//fprintf(stdout,"match score at position %i: %f\n",i,score);
		temp_best+=score;
		best_match_scores[i]=temp_best;
    
	}
  
	//fprintf(stdout,"Best match scores...END\n");

	/*********************************************************************************************/
	/*Left and right alignments */
	/*********************************************************************************************/
//	fprintf(stdout,"Left and right alignments...\n");

	std::vector<int> comp_sites;
	fast_fill_side_unspliced_first(nr_paths_par,seed_matrix_right,read_len,dna_len, read, dna, prb,functions, matchmatrix,qualityScores, donor,acceptor,comp_sites,seed_i, 
								   seed_j,best_match_scores,true,true,max_number_introns,max_gap,max_mism,max_edit_op,min_match, verbosity,currentMode,remapping);
//	fprintf(stdout,"%i right sides of the matrix filled...\n",seed_matrix_right.size());
  
	comp_sites.clear();
	fast_fill_side_unspliced_first(nr_paths_par,seed_matrix_left,read_len,dna_len, read, dna, prb,functions, matchmatrix,qualityScores, acceptor,donor,comp_sites,seed_i, 
								   seed_j,best_match_scores,false,true,max_number_introns,max_gap,max_mism,max_edit_op,min_match, verbosity,currentMode,remapping);
	comp_sites.clear();
//	fprintf(stdout,"%i left sides of the matrix filled...\n",seed_matrix_left.size());

//	fprintf(stdout,"Left and right alignments...END\n");


//   for(int n=0;n<seed_matrix_left.size();n++){
//     if (((SeedElem*)seed_matrix_left[n])!=NULL)
//       fprintf(stdout,"seed position %i %i %f\n",((SeedElem*)seed_matrix_left[n])->read_pos,((SeedElem*)seed_matrix_left[n])->dna_pos,((SeedElem*)seed_matrix_left[n])->best_scores[0]);
//   }


	/*********************************************************************************************/
	/* Find out the nr_paths_par best combinations */
	/*********************************************************************************************/

	memset(max_score_positions,0,sizeof(int)*2*nr_paths_par);
	if(seed_matrix_left[0]!=NULL && seed_matrix_right[0]!=NULL){ 
   
		int best_left_left=1;
		int best_left_right=0;
		int best_right_left=0;
		int best_right_right=1;

		max_score_positions[0]=0;
		max_score_positions[1]=0;

		for(int z=1;z<nr_paths_par;z++){
			double temp_left= seed_matrix_left[0]->best_scores[best_left_left] + seed_matrix_right[0]->best_scores[best_left_right];
			double temp_right= seed_matrix_left[0]->best_scores[best_right_left] + seed_matrix_right[0]->best_scores[best_right_right];

			if (temp_left > temp_right){
				max_score_positions[2*z]=best_left_left;
				max_score_positions[2*z+1]=best_left_right;
				best_left_left++;
			}
			else{
				max_score_positions[2*z]=best_right_left;
				max_score_positions[2*z+1]=best_right_right;
				best_right_right++;
			}
		}
	}

	/*********************************************************************************************/
	/* Display results */
	/*********************************************************************************************/

	// for (int z=0;z<nr_paths_par;z++){

	//   int z_path_left=max_score_positions[2*z] ; //path number for the left seed matrix
	//   int z_path_right=max_score_positions[2*z+1] ; //path number for the right seed matrix

	//   fprintf(stdout,"Align result for %i left matrix and %i right matrix \n",z_path_left,z_path_right);
  
	//   if(seed_matrix_left[0]!=NULL && seed_matrix_right[0]!=NULL && 
	//      seed_matrix_left[0]->best_scores[z_path_left]>-ALMOST_INFINITY && seed_matrix_right[0]->best_scores[z_path_right]>-ALMOST_INFINITY){

	//     fprintf(stdout,"score: %f\n",seed_matrix_left[0]->best_scores[z_path_left] + seed_matrix_right[0]->best_scores[z_path_right]);
	//     fprintf(stdout, "seed score %f\n",getScore(qualityScores,MMATRIX_LEN,check_char(read[seed_i]),check_char(dna[seed_j]),read_scores[seed_i]));

	//     SeedElem *next_seed=seed_matrix_left[0];
	//     int zz=z_path_left;

	//     int rstart;
	//     int dstart;
    
	//     std::vector<char> read_align;
	//     std::vector<char> dna_align;
	//     int max_gaps2;

	//     while (next_seed!=NULL){
	// 	std::vector<char> read_align_temp;
	// 	std::vector<char> dna_align_temp;      
	// 	max_gaps2=next_seed->max_gaps;
	// 	Prev_score** matrix=next_seed->matrices;
	// 	rstart=next_seed->best_score_pos[zz]->read_pos;
	// 	dstart=next_seed->best_score_pos[zz]->dna_pos;
	// 	int rseed=next_seed->read_pos;
	// 	int dseed=next_seed->dna_pos;
	// 	int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
      
	// 	while(!(rstart==rseed+1 && dstart==dseed+1)){
       
	// 	  //fprintf(stdout,"%i-%i:%i-%i (%i)\n",rstart,dstart,rseed,dseed,prev_z);
	// 	  int matrix_position= (rseed-rstart)*(max_gaps2*2+1)+(dseed-(rseed-rstart))-dstart+max_gaps2; 
	// 	  int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
	// 	  int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
	// 	  prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
	// 	  //fprintf(stdout,"prev %i-%i:%i-%i (%i)\n",rstart,dstart,prev_i,prev_j, prev_z);
	// 	  assert(rstart<=prev_i && dstart<=prev_j);
	  
	// 	  if (prev_i==rstart && prev_j==dstart+1){//read gap
	// 	    dna_align_temp.push_back(dna[dstart]);
	// 	    read_align_temp.push_back('-');	    
	// 	  }


	// 	  else if(prev_i==rstart+1 && prev_j==dstart){//dna gap
	// 	    dna_align_temp.push_back('-');
	// 	    read_align_temp.push_back(read[rstart]);
	// 	  }


	// 	  else if(prev_i==rstart+1 && prev_j==dstart+1){//match/mismatch
	// 	    dna_align_temp.push_back(dna[dstart]);
	// 	    read_align_temp.push_back(read[rstart]);
	// 	  }

	// 	  rstart=prev_i;
	// 	  dstart=prev_j;
	
	// 	}
	 
     
	// 	rstart=next_seed->best_score_pos[zz]->read_pos;
	// 	dstart=next_seed->best_score_pos[zz]->dna_pos;
	
	// 	int tmp=next_seed->best_score_pos[zz]->path_number;
	//     	next_seed=next_seed->best_score_pos[zz]->next_seed;
	// 	zz=tmp;
	
	// 	std::vector<char>::reverse_iterator rit;
	// 	for ( rit= dna_align_temp.rbegin() ; rit < dna_align_temp.rend(); ++rit )
	// 	  dna_align.push_back(*rit);
	// 	for ( rit= read_align_temp.rbegin() ; rit < read_align_temp.rend(); ++rit )
	// 	  read_align.push_back(*rit);
	
	// 	dna_align_temp.clear();
	// 	read_align_temp.clear();
	
	// 	if (next_seed!=NULL){
	// 	  for(int n=dstart-1;n>=next_seed->dna_pos+1;n--){
	// 	    dna_align.push_back(dna[n]);
	// 	    read_align.push_back('*');
	// 	  } 	   
	// 	}
	//     }
    
      
	//     std::vector<char>::reverse_iterator rit;
	//     fprintf(stdout,"DNA  ALIGN: ");
	//     for ( rit= dna_align.rbegin() ; rit < dna_align.rend(); ++rit )
	// 	fprintf(stdout,"%c",*rit);
	//     fprintf(stdout,"\nREAD ALIGN: ");
	//     for ( rit= read_align.rbegin() ; rit < read_align.rend(); ++rit )
	// 	fprintf(stdout,"%c",*rit);
	//     fprintf(stdout,"\n");
      
	//     dna_align.clear();
	//     read_align.clear();
      
	//     next_seed=seed_matrix_right[0];
	//     zz=z_path_right;

	//     while (next_seed!=NULL){
	// 	std::vector<char> read_align_temp;
	// 	std::vector<char> dna_align_temp;      
	// 	max_gaps2=next_seed->max_gaps;      
	// 	Prev_score** matrix=next_seed->matrices;
	// 	rstart=next_seed->best_score_pos[zz]->read_pos;
	// 	dstart=next_seed->best_score_pos[zz]->dna_pos;
	// 	int rseed=next_seed->read_pos;
	// 	int dseed=next_seed->dna_pos;
	// 	int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
      
	// 	while(!(rstart==rseed-1 && dstart==dseed-1)){
	
	// 	  int matrix_position= (rstart-rseed)*(max_gaps2*2+1)+dstart-(dseed+rstart-rseed)+max_gaps2;
	// 	  //	fprintf(stdout,"(%i-%i),(%i,%i) %i\n",rstart,dstart,rseed,dseed,matrix_position);
	// 	  int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
	// 	  int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
	// 	  prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
	// 	  //fprintf(stdout,"prev %i-%i:%i-%i (%i)\n",rstart,dstart,prev_i,prev_j, prev_z);
	  
	// 	  if (prev_i==rstart && prev_j==dstart-1){//read gap
	// 	    dna_align_temp.push_back(dna[dstart]);
	// 	    read_align_temp.push_back('-');
	// 	  }
	// 	  else if(prev_i==rstart-1 && prev_j==dstart){//dna gap
	// 	    dna_align_temp.push_back('-');
	// 	    read_align_temp.push_back(read[rstart]);
	// 	  }
	// 	  else if(prev_i==rstart-1 && prev_j==dstart-1){//match/mismatch
	// 	    dna_align_temp.push_back(dna[dstart]);
	// 	    read_align_temp.push_back(read[rstart]);
	// 	  }

	
	// 	  rstart=prev_i;
	// 	  dstart=prev_j;
	// 	}
      
	
	// 	rstart=next_seed->best_score_pos[zz]->read_pos;
	// 	dstart=next_seed->best_score_pos[zz]->dna_pos;

	// 	int tmp=next_seed->best_score_pos[zz]->path_number;
	// 	next_seed=next_seed->best_score_pos[zz]->next_seed;
	// 	zz=tmp;


	// 	for ( rit= dna_align_temp.rbegin() ; rit < dna_align_temp.rend(); ++rit )
	// 	  dna_align.push_back(*rit);
	// 	for ( rit= read_align_temp.rbegin() ; rit < read_align_temp.rend(); ++rit )
	// 	  read_align.push_back(*rit);
	
	// 	dna_align_temp.clear();
	// 	read_align_temp.clear();
	      
	// 	if (next_seed!=NULL){
	// 	  for(int n=dstart+1;n<=next_seed->dna_pos-1;n++){
	// 	    dna_align.push_back(dna[n]);
	// 	    read_align.push_back('*');
	// 	  } 
	// 	}
	//     }


	//     std::vector<char>::iterator it;
	//     fprintf(stdout,"DNA  ALIGN: ");
	//     for ( it= dna_align.begin() ; it < dna_align.end(); it++ )
	// 	fprintf(stdout,"%c",*it);
	//     fprintf(stdout,"\nREAD ALIGN: ");
	//     for ( it= read_align.begin() ; it < read_align.end(); it++ )
	// 	fprintf(stdout,"%c",*it);
	//     fprintf(stdout,"\n");
	//     dna_align.clear();
	//     read_align.clear();
	//   }
	// }
  
	// SeedElem* next_seed= seed_matrix_left[0];
	// while(next_seed!=NULL){
	//   fprintf(stdout,"seed %i-%i\n",next_seed->read_pos,next_seed->dna_pos);
	//   next_seed=next_seed->best_score_pos[0]->next_seed;
	// }


	/*********************************************************************************************/
	/* Clean structures */
	/*********************************************************************************************/
	//  fprintf(stdout,"Clean structures...\n");
 
	delete[] best_match_scores;
	best_match_scores=NULL;

	//  fprintf(stdout,"Clean structures...END\n");
}



