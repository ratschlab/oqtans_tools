// Authors: Bettina Hepp, Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean, Soeren Sonnenburg
// Copyright (C) 2005-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include "qpalma_dp.h"

#include <cstring>

using namespace std;

/*[splice_align, est_align, weightMatch, alignmentscores, dnaest] = ...
  myalign([nr_paths], dna, est, {h}, matchmatrix, donor, acceptor, remove_duplicate_scores, ...
  min_match) */
Alignment::Alignment(int numQPlifs, int numq, bool use_qscores, int _verbosity) {
      len = 0;
      limits = 0;
      penalties = 0;
      max_len = 0;
      min_len = 0;
      cache = 0;
      transform = T_LINEAR;
      id = 0;
      name = 0;
      use_svm = 0;
	  DNA_ARRAY = NULL ;
	  EST_ARRAY = NULL ;
	  verbosity=_verbosity ;

      // set ptrs to zero first
      splice_align         = 0;
      est_align            = 0;
      mmatrix_param        = 0;
      alignmentscores      = 0;
      qualityFeaturesAllPaths = 0;
      mlen = 6; // score matrix: length of 6 for "- A C G T N"

      numQualSuppPoints = numq;
      numPlifs = numQPlifs;
      use_quality_scores = use_qscores;

      //printf("number of support points: %d\n",numQualSuppPoints);
      //printf("number of plifs: %d\n",numPlifs );
      assert( numQualSuppPoints >= 0 );
      assert( numPlifs >= 0 );
}

void Alignment::getDNAEST(){}

/* 
   est -> read sequence
   prb -> read quality values
   h -> intron length scoring function
   matchmatrix -> 
   mm_len -> size of matchmatrix
   donor -> scores for donor splice site prediction for every position in dna (-inf if not donor splice site)
   d_len -> length(donor)
   acceptor -> scores  for acceptor splice site prediction for every position in dna (-inf if not acceptor splice site)
   qualityScores -> 
   remove_duplicate_scores -> hack -> false
*/
//Add starting position for alignment according to a long hit (read and dna pos, len) and strand
void Alignment::myalign_fast(char strand, Chromosome const &chr,  std::vector<int> &positions, int nr_paths_p, char* dna, int dna_len_p, char* est,    int est_len_p, double* prb, struct penalty_struct h, double* matchmatrix, int mm_len,
							 double* donor, int d_len, double* acceptor, int a_len, struct penalty_struct* qualityScores, 
							 bool remove_duplicate_scores, int hit_read, int hit_dna, int hit_len, int max_number_introns,  
							 int max_gap, int max_mism, int max_edit_op, int min_match, bool remapping) {

	// printf("Entering myalign_fast...\n");
	nr_paths = nr_paths_p;
	dna_len=dna_len_p;
	est_len=est_len_p;

	mode currentMode;
	if (use_quality_scores) 
		currentMode = USE_QUALITY_SCORES;
	else
		currentMode = NORMAL;

	// dnaest
	DNA_ARRAY = 0;
	EST_ARRAY = 0;


	/***************************************************************************/ 
	// initialize seed_position, seed vectors and call fast_fill_matrix()  
	/***************************************************************************/

	///////////////////////////////
	// Initialize seed positions for left and right alignments
	///////////////////////////////
	int seed_i=hit_read; 
	int seed_j=hit_dna; 

	double best_match=-ALMOST_INFINITY;
	const int end_error=5;

	for(int i=0+end_error; i < hit_len-end_error && hit_read+i<est_len_p && hit_dna+i<dna_len_p ;i++){
    
		int i_pos=hit_read+i ;
		int j_pos=hit_dna+i ;
		if (i_pos<0 || i_pos>=est_len_p)
			continue ;
		if (j_pos<0 || j_pos>=dna_len_p)
			continue ;
	
		double score;
		if (currentMode == USE_QUALITY_SCORES)
			score=getScore(qualityScores,mm_len,check_char(dna[j_pos]),check_char(est[i_pos]),prb[i_pos]) ;
		else
			score=matchmatrix[mm_len*check_char(dna[j_pos])+check_char(est[i_pos])] ;

		if (score>best_match){
			best_match=score;
			seed_i=i_pos;
			seed_j=j_pos;
		}
	}
  
	// fprintf(stdout,"seed (%i-%i) %f\n",seed_i,seed_j,best_match);



	///////////////////////////////
	// Initialize seed vectors for left and right alignments
	///////////////////////////////
	std::vector<SeedElem *> seed_matrix_left;
	std::vector<SeedElem *> seed_matrix_right;
	int* max_score_positions = new int[nr_paths*2];

	fast_fill_matrix(nr_paths, max_score_positions, est_len, dna_len, est, dna, prb, &h, matchmatrix, qualityScores, 
					 donor, acceptor,remove_duplicate_scores,seed_i,seed_j,seed_matrix_left, seed_matrix_right, max_number_introns,max_gap,max_mism,max_edit_op,min_match, verbosity,currentMode,remapping);


	/***************************************************************************/ 
	// return arguments etc.  
	/***************************************************************************/ 
	int result_length; 
  
	splice_align_size = (dna_len)*nr_paths;
	est_align_size = (est_len)*nr_paths;

	int mmatrix_size;

	if (currentMode == USE_QUALITY_SCORES) {
		mmatrix_param_size = mlen*nr_paths;
		mmatrix_size = mlen;
	}
  
	if (currentMode == NORMAL) {
		mmatrix_param_size = (mlen*mlen)*nr_paths;
		mmatrix_size = mlen*mlen;
	}
   
	alignmentscores_size = nr_paths; //alignment score for each path/matrix
	numPathsPlifs = numPlifs*nr_paths; //alignment score for each path/matrix

	splice_align         = new int[splice_align_size];
	est_align            = new int[est_align_size];
	mmatrix_param        = new int[mmatrix_param_size];
	alignmentscores      = new double[nr_paths]; //alignment score for each path/matrix


	//  printf("before memset...\n");
	memset((char*)splice_align, -1, (dna_len)*nr_paths*sizeof(int)); // fills splice_align with zeros
	memset((char*)est_align, -1, (est_len)*nr_paths*sizeof(int)); // fills est_align with zeros
	memset((char*)mmatrix_param, 0, mmatrix_size*nr_paths*sizeof(int)); //fills mmatrix_param with zeros
	memset(alignmentscores, -1, nr_paths*sizeof(double)); //fills alignmentscores with zeros
	//printf("after memset...\n");
  
	qualityFeaturesAllPaths= new penalty_struct*[nr_paths];

	for (int z=0; z<nr_paths; z++) {
		result_length = 0 ;
    
		int* s_align = splice_align + (dna_len)*z;  //pointer
		int* e_align = est_align + (est_len)*z ;    //pointer
		int* mparam  = mmatrix_param + mmatrix_size*z; //pointer

		qualityFeaturesAllPaths[z] = new penalty_struct[numPlifs];
    
		int qidx, pidx;
		for(qidx=0;qidx<numPlifs;qidx++) {
			penalty_struct p;
			init_penalty_struct(p);
			p.len = numQualSuppPoints;
			p.limits = (double*) calloc(p.len,sizeof(double));

			for(pidx=0;pidx<p.len;pidx++)
				p.limits[pidx] = qualityScores[qidx].limits[pidx];

			p.penalties = (double*) calloc(p.len,sizeof(double));
			qualityFeaturesAllPaths[z][qidx] = p;
		}
	
		//penalty_struct* qparam = qualityScoresAllPaths + (numPlifs*z);

		//    printf("before call to fast_result_align...\n");
		bool no_more_path = fast_result_align(seed_matrix_left,seed_matrix_right, z, est_len, dna_len, &result_length, est, dna, prb, s_align, e_align, mparam, alignmentscores, max_score_positions, qualityFeaturesAllPaths[z] , currentMode,best_match);
		//printf("after call to fast_result_align...\n");
	
		//printf("z is %d\n",z);
		//int len;
		//for(qidx=0;qidx<numPlifs;qidx++) {
		//   penalty_struct p;
		//   p = qualityScoresAllPaths[z][qidx];
		//   printf("%d: ",qidx);
		//   for(pidx=0;pidx<p.len;pidx++)
		//      printf("%f ",p.limits[pidx]);
		//   printf("\n");

		//   for(pidx=0;pidx<p.len;pidx++)
		//      printf("%f ",p.penalties[pidx]);
		//   printf("\n");
		//}

		if(z==0) {

			if(DNA_ARRAY != 0) {
				delete[] DNA_ARRAY;
				delete[] EST_ARRAY;
			}
		
			result_len = result_length;
     
			DNA_ARRAY = new int[result_length];
			EST_ARRAY = new int[result_length];

			if (no_more_path==0){

				int z_path_left=max_score_positions[2*z] ; //path number for the left seed matrix
				int z_path_right=max_score_positions[2*z+1] ; //path number for the right seed matrix
     
         
       
				SeedElem *next_seed=seed_matrix_left[0];
				int zz=z_path_left;
     
				std::vector<char> read_align;
				std::vector<char> dna_align;
       
				while (next_seed!=NULL){
	 
      
					std::vector<char> read_align_temp;
					std::vector<char> dna_align_temp;
	 
					Prev_score** matrix=next_seed->matrices;
					int rstart=next_seed->best_score_pos[zz]->read_pos;
					int dstart=next_seed->best_score_pos[zz]->dna_pos;
					int rseed=next_seed->read_pos;
					int dseed=next_seed->dna_pos;
					int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
					int gaps_for_seed=next_seed->max_gaps;

					while(!(rstart==rseed+1 && dstart==dseed+1)){
	 
						int matrix_position= (rseed-rstart)*(gaps_for_seed*2+1)+(dseed-(rseed-rstart))-dstart+gaps_for_seed;
						int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
						int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
						prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
						//fprintf(stdout,"prev %i-%i:%i-%i (%i)\n",rstart,dstart,prev_i,prev_j, prev_z);

						assert(rstart<=prev_i && dstart<=prev_j);
	  
						if (prev_i==rstart && prev_j==dstart+1){//read gap
							if (dna[dstart]=='N'){
								if(strand=='-')
									dna_align_temp.push_back(check_char(get_compl_base(chr[positions[dstart]])));
								else
									dna_align_temp.push_back(check_char(chr[positions[dstart]]));
							}
							else
								dna_align_temp.push_back(check_char(dna[dstart]));
							read_align_temp.push_back(0);	    
						}

						else if(prev_i==rstart+1 && prev_j==dstart){//dna gap
							dna_align_temp.push_back(0);
							read_align_temp.push_back(check_char(est[rstart]));	    
						}


						else if(prev_i==rstart+1 && prev_j==dstart+1){//match/mismatch
							if (dna[dstart]=='N'){
								if(strand=='-')
									dna_align_temp.push_back(check_char(get_compl_base(chr[positions[dstart]])));
								else
									dna_align_temp.push_back(check_char(chr[positions[dstart]]));
							}
							else
								dna_align_temp.push_back(check_char(dna[dstart]));
							read_align_temp.push_back(check_char(est[rstart]));	    
						}
	   
						rstart=prev_i;
						dstart=prev_j;
	   
					}
	 
					rstart=next_seed->best_score_pos[zz]->read_pos;
					dstart=next_seed->best_score_pos[zz]->dna_pos;
	 
					int tmp=next_seed->best_score_pos[zz]->path_number;
					next_seed=next_seed->best_score_pos[zz]->next_seed;
					zz=tmp;
	 
					std::vector<char>::reverse_iterator rit;
					for ( rit= dna_align_temp.rbegin() ; rit < dna_align_temp.rend(); ++rit )
						dna_align.push_back(*rit);
					for ( rit= read_align_temp.rbegin() ; rit < read_align_temp.rend(); ++rit )
						read_align.push_back(*rit);
	 
					dna_align_temp.clear();
					read_align_temp.clear();
	 
					if (next_seed!=NULL)
						for(int n=dstart-1;n>=next_seed->dna_pos+1;n--){
							dna_align.push_back(check_char(dna[n]));
							read_align.push_back(6);
						} 
				}
     
     
				int size=dna_align.size();
				for(int n=0; n<size;n++){
       
					DNA_ARRAY[size-n-1]=dna_align[n];
					EST_ARRAY[size-n-1]=read_align[n];
				}
     
				dna_align.clear();
				read_align.clear();
     
				next_seed=seed_matrix_right[0];
				zz=z_path_right;

     
				while (next_seed!=NULL){

					std::vector<char> read_align_temp;
					std::vector<char> dna_align_temp;      
					Prev_score** matrix=next_seed->matrices;
					int rstart=next_seed->best_score_pos[zz]->read_pos;
					int dstart=next_seed->best_score_pos[zz]->dna_pos;
					int rseed=next_seed->read_pos;
					int dseed=next_seed->dna_pos;
					int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
					int gaps_for_seed=next_seed->max_gaps;
       
					while(!(rstart==rseed-1 && dstart==dseed-1)){
	
						int matrix_position= (rstart-rseed)*(2*gaps_for_seed+1)+dstart-(dseed+rstart-rseed)+gaps_for_seed;
						//	fprintf(stdout,"(%i-%i),(%i,%i) %i\n",rstart,dstart,rseed,dseed,matrix_position);
						int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
						int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
						prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
	  
						if (prev_i==rstart && prev_j==dstart-1){//read gap
							if (dna[dstart]=='N'){
								if(strand=='-')
									dna_align_temp.push_back(check_char(get_compl_base(chr[positions[dstart]])));
								else
									dna_align_temp.push_back(check_char(chr[positions[dstart]]));
							}
							else
								dna_align_temp.push_back(check_char(dna[dstart]));
							read_align_temp.push_back(0);
						}
						else if(prev_i==rstart-1 && prev_j==dstart){//dna gap
							dna_align_temp.push_back(0);
							read_align_temp.push_back(check_char(est[rstart]));
						}
						else if(prev_i==rstart-1 && prev_j==dstart-1){//match/mismatch
							if (dna[dstart]=='N'){
								if(strand=='-')
									dna_align_temp.push_back(check_char(get_compl_base(chr[positions[dstart]])));
								else
									dna_align_temp.push_back(check_char(chr[positions[dstart]]));
							}
							else
								dna_align_temp.push_back(check_char(dna[dstart]));
							read_align_temp.push_back(check_char(est[rstart]));
						}

						rstart=prev_i;
						dstart=prev_j;
	
					}
      
	
					rstart=next_seed->best_score_pos[zz]->read_pos;
					dstart=next_seed->best_score_pos[zz]->dna_pos;
       
					int tmp=next_seed->best_score_pos[zz]->path_number;
					next_seed=next_seed->best_score_pos[zz]->next_seed;
					zz=tmp;

					std::vector<char>::reverse_iterator rit;
					for ( rit= dna_align_temp.rbegin() ; rit < dna_align_temp.rend(); ++rit )
						dna_align.push_back(*rit);
					for ( rit= read_align_temp.rbegin() ; rit < read_align_temp.rend(); ++rit )
						read_align.push_back(*rit);
       
					dna_align_temp.clear();
					read_align_temp.clear();
	
					if (next_seed!=NULL){
	   
						for (int n=dstart+1;n<=next_seed->dna_pos-1;n++){
							dna_align.push_back(check_char(dna[n]));
							read_align.push_back(6);
						} 
					} 
				}
  
				int size2=dna_align.size();
				for(int n=0; n<size2;n++){
	 
					DNA_ARRAY[(result_length-size2)+n]=dna_align[n];
					EST_ARRAY[(result_length-size2)+n]=read_align[n];
				}
       
				dna_align.clear();
				read_align.clear();

				// for (int n=0;n<result_len;n++)
				// 	fprintf(stdout,"%i",DNA_ARRAY[n]);
				// fprintf(stdout,"\n");
				// for (int n=0;n<result_len;n++)
				// 	fprintf(stdout,"%i",EST_ARRAY[n]);
				// fprintf(stdout,"\n");

			}
		} // end of "if z == 0"


	} //end of z



	/***************************************************************************/ 
	// Clean structures
	/***************************************************************************/
	clean_seed_matrix_vector(seed_matrix_right,nr_paths);
	clean_seed_matrix_vector(seed_matrix_left,nr_paths);
	delete[] max_score_positions ;

}


void Alignment::cleanup() 
{
	if (qualityFeaturesAllPaths)
	{
		for (int z=0; z<nr_paths; z++) 
		{
			for(int qidx=0;qidx<numPlifs;qidx++)
			{
				free(qualityFeaturesAllPaths[z][qidx].limits) ;
				free(qualityFeaturesAllPaths[z][qidx].penalties) ;
				delete[] qualityFeaturesAllPaths[z][qidx].cache ;
			}
			delete[] qualityFeaturesAllPaths[z] ;
		}
		
		delete[] qualityFeaturesAllPaths ;
		qualityFeaturesAllPaths=NULL ;
	}
	if(splice_align != NULL)
	{
		delete[] splice_align;
		splice_align=NULL ;
	}

	if(est_align != NULL)
	{
		delete[] est_align;
		est_align = NULL ;
	}
	
	if(mmatrix_param != NULL)
	{
		delete[] mmatrix_param;
		mmatrix_param=NULL ;
	}

	if(alignmentscores != NULL)
	{
		delete[] alignmentscores;
		alignmentscores=NULL ;
	}	
	
	if (DNA_ARRAY != NULL) {
		delete[] DNA_ARRAY;
		delete[] EST_ARRAY;
    }

}

void Alignment::getAlignmentResults(int* s_align, int* e_align,
      int* mmatrix_p, double* alignscores, double* qScores) {

   int idx;
   for(idx=0; idx<(int)splice_align_size; idx++){
      s_align[idx] = splice_align[idx];
   }
   for(idx=0; idx<(int)est_align_size; idx++)
      e_align[idx] =  est_align[idx];

   for(idx=0; idx<(int)mmatrix_param_size; idx++)
      mmatrix_p[idx] = mmatrix_param[idx];

   for(idx=0; idx<(int)alignmentscores_size; idx++)
      alignscores[idx] = alignmentscores[idx];
   
   if (use_quality_scores && qScores!=NULL) {
	   penalty_struct currentPlif;
	   int ctr=0;
	   for (int z=0; z<nr_paths; z++) {
	     for(int estChar=1;estChar<6;estChar++) {
	       for(int dnaChar=0;dnaChar<6;dnaChar++) {
		 
		 int currentPos = (estChar-1)*6+dnaChar;
		 currentPlif = qualityFeaturesAllPaths[z][currentPos];
		 
		 for(int pidx=0; pidx<currentPlif.len; pidx++) {
		   qScores[ctr] = currentPlif.penalties[pidx];
		   //printf("%f ",qScores[ctr]);
		   ctr++;
		 }
		 //printf("\n");
	       }}}
	   
	   //printf("\nctr is %d\n",ctr);
   }
   //printf("Leaving getAlignmentResults...\n");
}

double Alignment::scoreUnsplicedAlignment(const char * align_seq, double * prb, int read_length, struct penalty_struct* qualityScores, double * matchmatrix, char strand) 
{
	int len_=strlen(align_seq) ;
	double score=0.0 ;
	//fprintf(stdout, "align_seq=%s\n", align_seq) ;
	
	assert(strand=='+'||strand=='-') ;
	int reverse[7] = { -1, 0, 4, 3, 2, 1, 5 } ;
			
	int pos=0 ;
	for (int i=0; i<len_; i++)
	{
		int dnachar ;
		int estchar ;

		assert(align_seq[i]!=']') ;
		if (align_seq[i]=='[')
		{
			dnachar=check_char(align_seq[i+1]) ;
			estchar=check_char(align_seq[i+2]) ;
			assert(align_seq[i+3]==']') ;
			i+=3 ;
		}
		else
			estchar=dnachar=check_char(align_seq[i]) ;

		if (strand=='-')
		{
			assert(dnachar>=-1 && dnachar<=5) ;
			assert(estchar>=-1 && estchar<=5) ;
			
			dnachar = reverse[dnachar+1] ;
			estchar = reverse[estchar+1] ;
		}

		if (dnachar==-1)
			dnachar=0 ;

		//fprintf(stderr, "i=%i, pos=%i, dnachar=%i, estchar=%i\n", i, pos, dnachar, estchar) ;
		double score_ ;
		if (estchar>0)
		{
			assert(pos<read_length) ;
			if (use_quality_scores)
				score_=getScore(qualityScores, 6, dnachar, estchar, prb[pos]) ;
			else
				score_=matchmatrix[6*dnachar+estchar] ;
		}
		else
			score_=matchmatrix[dnachar] ;
				
		//fprintf(stderr, "looking up %i %i at %1.1f -> %1.2f", dnachar, estchar, prb[pos], score_) ;
		
		score+=score_ ;
		
		if (estchar>0)
			pos++ ;
	}
	//assert(pos==read_length);
	
	return score ;
}

void Alignment::getAlignmentArrays(int* dna_align, int* read_align) {

   int idx;
   for(idx=0; idx<result_len; idx++) {
      dna_align[idx] = DNA_ARRAY[idx];
      read_align[idx] = EST_ARRAY[idx];
   }
}
