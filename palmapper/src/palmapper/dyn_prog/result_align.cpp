#include "assert.h"
#include "fill_matrix.h"
#include "debug_tools.h"

using namespace std;

void increaseFeatureCount(penalty_struct* qparam, int dnaChar, int estChar, double estprb) {

   assert(estChar > 0 && estChar < 6);
   assert(dnaChar > -1 && dnaChar < 6);

   int currentPos = (estChar-1)*6+dnaChar;

   penalty_struct currentStruct = qparam[currentPos];

   //printf("Current index %d est/dna: %d %d score %f\n",currentPos,estChar,dnaChar,estprb);

   //if( estChar == 1 and dnaChar == 1) {
   //   printf("before\n");
   //   int p_idx;
   //   for(p_idx=0;p_idx<currentStruct.len;p_idx++) {
   //      printf("%f ",currentStruct.limits[p_idx]);
   //   }
   //   printf("\n");

   //   for(p_idx=0;p_idx<currentStruct.len;p_idx++) {
   //      printf("%f ",currentStruct.penalties[p_idx]);
   //   }
   //   printf("\n");
   //}

   double value = estprb;
   int Lower = 0;
   int idx;

   for (idx=0;idx<currentStruct.len;idx++) {
      if (currentStruct.limits[idx] <= value)
         Lower++;
   }

   if (Lower == 0) {
         currentStruct.penalties[0] += 1;
         qparam[currentPos] = currentStruct;
         return;
   }

   if (Lower == currentStruct.len) {
         currentStruct.penalties[currentStruct.len-1] += 1;
         qparam[currentPos] = currentStruct;
         return;
   }

   Lower -= 1;
   int Upper = Lower+1; // x-werte bleiben fest

   double weightup  = 1.0*(value - currentStruct.limits[Lower]) / (currentStruct.limits[Upper] - currentStruct.limits[Lower]);
   double weightlow = 1.0*(currentStruct.limits[Upper] - value) / (currentStruct.limits[Upper] - currentStruct.limits[Lower]);
   currentStruct.penalties[Upper] += weightup;
   currentStruct.penalties[Lower] += weightlow; 

   //if( estChar == 1 and dnaChar == 1) {
   //   printf("estprb/Lower/Upper %f %d %d\n",estprb,Lower,Upper);
   //   int p_idx;
   //   printf("after\n");
   //   for(p_idx=0;p_idx<currentStruct.len;p_idx++) {
   //      printf("%f ",currentStruct.limits[p_idx]);
   //   }
   //   printf("\n");
   //   for(p_idx=0;p_idx<currentStruct.len;p_idx++)
   //      printf("%f ",currentStruct.penalties[p_idx]);
   //   printf("\n");
   //}
   qparam[currentPos] = currentStruct;
}



bool fast_result_align(const std::vector<SeedElem*>& seed_matrix_left,const std::vector<SeedElem*>& seed_matrix_right, int z, int est_len, int dna_len, int* result_length_ptr, 
			char* est, char* dna, double* prb, int* s_align, int* e_align, int* mparam, double* alignmentscores, int* max_score_positions, 
			penalty_struct* qparam, mode currentMode, double score_seed){  

  const int mlen=6; // length of matchmatrix

  int z_path_left=max_score_positions[2*z] ; //path number for the left seed matrix
  int z_path_right=max_score_positions[2*z+1] ; //path number for the right seed matrix

  int dnanum ; //using check_char: letter -> number
  int estnum ; //using check_char: letter -> number

  double prbnum ;
  
  int splice_state ;
  int est_state ;

  splice_state=0;
  est_state=1;

  alignmentscores[z] = -ALMOST_INFINITY;

  //  fprintf(stdout,"Align result for path %i of left matrix and path %i of right matrix \n",z_path_left,z_path_right);
  
  //First elements in seed vectors have to be not null
  //Scores for right and left paths have to be not -inf
  if(seed_matrix_left[0]!=NULL && seed_matrix_right[0]!=NULL && 
     seed_matrix_left[0]->best_scores[z_path_left]>-ALMOST_INFINITY &&seed_matrix_right[0]->best_scores[z_path_right]>-ALMOST_INFINITY){

    // Alignment score: the score at the seed position was counted twice (once for each side)
    alignmentscores[z] = seed_matrix_left[0]->best_scores[z_path_left] + seed_matrix_right[0]->best_scores[z_path_right] - score_seed;


    /**********************************************************/
    /* LEFT SIDE                                              */
    /**********************************************************/
    //Init next_seed from the seed position
    SeedElem *next_seed=seed_matrix_left[0];
    int zz=z_path_left;

    int rstart=0;
    int dstart=0;
    int max_gap;

    while (next_seed!=NULL){
      
      //Number of allowed gaps in this matrix (determines the width of the band)
      max_gap=next_seed->max_gaps;
      Prev_score** matrix=next_seed->matrices;
      rstart=next_seed->best_score_pos[zz]->read_pos;
      dstart=next_seed->best_score_pos[zz]->dna_pos;
      int rseed=next_seed->read_pos;
      int dseed=next_seed->dna_pos;
      int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
      
      while(!(rstart==rseed+1 && dstart==dseed+1)){
		  
		  //fprintf(stdout,"%i-%i:%i-%i (%i)\n",rstart,dstart,rseed,dseed,prev_z);
		  int matrix_position= (rseed-rstart)*(max_gap*2+1)+(dseed-(rseed-rstart))-dstart+max_gap; 
		  int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
		  int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
		  prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
		  //	fprintf(stdout,"prev %i-%i:%i-%i (%i)\n",rstart,dstart,prev_i,prev_j, prev_z);
		  assert(rstart<=prev_i && dstart<=prev_j);
		  
		  //read gap
		  if (prev_i==rstart && prev_j==dstart+1){
			  (*result_length_ptr)= (*result_length_ptr) + 1;
			  s_align[dstart] = splice_state; 
			  
			  dnanum = check_char(dna[dstart]) ; 
			  estnum = 0 ; //gap
			  
			  if(currentMode == USE_QUALITY_SCORES)
				  mparam[dnanum] ++ ;
			  else
				  mparam[mlen*dnanum +estnum] ++ ;
		  }
		  
		  
		  //dna gap
		  else if(prev_i==rstart+1 && prev_j==dstart){
			  (*result_length_ptr)= (*result_length_ptr) + 1;
			  e_align[rstart] = est_state ;  //1 or 2, depended

			  dnanum = 0 ; //gap
			  estnum = check_char(est[rstart]) ;
			  
			  if(currentMode == USE_QUALITY_SCORES){
				  prbnum = prb[rstart];
				  increaseFeatureCount(qparam,dnanum,estnum,prbnum);
			  }
			  else
				  mparam[mlen*dnanum +estnum] ++;
		  }
		  
		  //match/mismatch
		  else if(prev_i==rstart+1 && prev_j==dstart+1){
			  (*result_length_ptr)= (*result_length_ptr) + 1;
			  
			  s_align[dstart] = splice_state; 
			  e_align[rstart] = est_state ; //1 or 2, depended
			  
			  dnanum = check_char(dna[dstart]); 
			  estnum = check_char(est[rstart]); 
			  
			  
			  if(currentMode == USE_QUALITY_SCORES){
				  prbnum = prb[rstart];
				  increaseFeatureCount(qparam,dnanum,estnum,prbnum);
			  }
			  else
				  mparam[mlen*dnanum+estnum] ++ ;
		  }

		  rstart=prev_i;
		  dstart=prev_j;
		  
      }
	      
      rstart=next_seed->best_score_pos[zz]->read_pos;
      dstart=next_seed->best_score_pos[zz]->dna_pos;
	  
      int tmp=next_seed->best_score_pos[zz]->path_number;
      next_seed=next_seed->best_score_pos[zz]->next_seed;
      zz=tmp;

      
      //Spliced alignment
      if (next_seed!=NULL){
		  //	  fprintf(stdout,"%i %i\n",dstart,next_seed->dna_pos-1);
		  (*result_length_ptr) =  (*result_length_ptr) + (dstart-next_seed->dna_pos-1);
		  
		  if (est_state==1) //was exon labeled "1"
			  est_state = 2; //new exon is labeled "2"
		  else
			  est_state = 1 ; //last exon was labeled "2", new exon is labeled "1"
		  
		  for (int n=dstart-1;n>=next_seed->dna_pos+1;n--){
			  if (splice_state == 0) //coming from exon
				  splice_state = 2; //first intron_state for left side: acceptor
			  else
				  splice_state = 3; //intron
			  
			  if (n == next_seed->dna_pos+1) //last intron_state for left side: donor
	    splice_state = 1;//donor
			  
			  s_align[n] = splice_state; 
		  }
		  
		  splice_state = 0 ; //exon again
      }
    }

    //Left side not aligned
    for (int dna_pos=0; dna_pos<dstart; dna_pos++) {
		s_align[dna_pos] = 4; 
    }
    for (int est_pos=0; est_pos<rstart; est_pos++) {
		e_align[est_pos] = 4;
    }            
    
    (*result_length_ptr) =  (*result_length_ptr) -1; //Seed position counted twice

    

    /**********************************************************/
    /* RIGHT SIDE                                             */
    /**********************************************************/
    //Init next_seed from the seed position
    next_seed=seed_matrix_right[0];
    zz=z_path_right;
    splice_state=0;
    est_state=1;

    while (next_seed!=NULL){
      
      max_gap=next_seed->max_gaps;      
      Prev_score** matrix=next_seed->matrices;
      rstart=next_seed->best_score_pos[zz]->read_pos;
      dstart=next_seed->best_score_pos[zz]->dna_pos;
      int rseed=next_seed->read_pos;
      int dseed=next_seed->dna_pos;
      int prev_z= next_seed->best_score_pos[zz]->path_number_matrices;
      
      while(!(rstart==rseed-1 && dstart==dseed-1)){
	
		  int matrix_position= (rstart-rseed)*(max_gap*2+1)+dstart-(dseed+rstart-rseed)+max_gap;
		  //fprintf(stdout,"(%i-%i),(%i,%i) %i\n",rstart,dstart,rseed,dseed,matrix_position);
		  int prev_i=((Prev_score*)matrix[prev_z]+matrix_position)->prev_i;
		  int prev_j=((Prev_score*)matrix[prev_z]+matrix_position)->prev_j;
		  prev_z=((Prev_score*)matrix[prev_z]+matrix_position)->prev_matrix_no;
		  //fprintf(stdout,"prev %i-%i:%i-%i (%i)\n",rstart,dstart,prev_i,prev_j, prev_z);
	  
		  //read gap
		  if (prev_i==rstart && prev_j==dstart-1){
			  (*result_length_ptr)= (*result_length_ptr) + 1;
			  s_align[dstart] = splice_state; 
	  
			  dnanum = check_char(dna[dstart]) ; 
			  estnum = 0 ; //gap

			  if(currentMode == USE_QUALITY_SCORES)
				  mparam[dnanum] ++ ;
			  else
				  mparam[mlen*dnanum +estnum] ++ ;
	  
		  }

		  //dna gap
		  else if(prev_i==rstart-1 && prev_j==dstart){
			  (*result_length_ptr)= (*result_length_ptr) + 1;
			  e_align[rstart] = est_state ;  //1 or 2, depended

			  dnanum = 0 ; //gap
			  estnum = check_char(est[rstart]) ;

			  if(currentMode == USE_QUALITY_SCORES){
				  prbnum = prb[rstart];
				  increaseFeatureCount(qparam,dnanum,estnum,prbnum);
			  }
			  else
				  mparam[mlen*dnanum +estnum] ++;
		  }

		  //match/mismatch
		  else if(prev_i==rstart-1 && prev_j==dstart-1){
			  (*result_length_ptr)= (*result_length_ptr) + 1;
      
			  s_align[dstart] = splice_state; 
			  e_align[rstart] = est_state ; //1 or 2, depended
	  
			  dnanum = check_char(dna[dstart]); 
			  estnum = check_char(est[rstart]); 
	  
			  if(currentMode == USE_QUALITY_SCORES){
				  prbnum = prb[rstart];
				  increaseFeatureCount(qparam,dnanum,estnum,prbnum);
			  }
			  else
				  mparam[mlen*dnanum+estnum] ++ ;
		  }

	
		  rstart=prev_i;
		  dstart=prev_j;
      }
      
	
      rstart=next_seed->best_score_pos[zz]->read_pos;
      dstart=next_seed->best_score_pos[zz]->dna_pos;

      int tmp=next_seed->best_score_pos[zz]->path_number;
      next_seed=next_seed->best_score_pos[zz]->next_seed;
      zz=tmp;

      //Spliced alignment
      if (next_seed!=NULL){
		  (*result_length_ptr) =  (*result_length_ptr) + (next_seed->dna_pos-dstart-1);
		  if (est_state==1) //was exon labeled "1"
			  est_state = 2; //new exon is labeled "2"
		  else
			  est_state = 1 ; //last exon was labeled "2", new exon is labeled "1"
	  
		  for (int n=dstart+1;n<=next_seed->dna_pos-1;n++){
			  if (splice_state == 0) //coming from exon
				  splice_state = 1; //first intron_state for right side: donor
			  else
				  splice_state = 3; //intron
	    
			  if (n == next_seed->dna_pos-1) //last intron_state for right side: acceptor
				  splice_state = 1;//acceptor
	  
			  s_align[n] = splice_state; 
		  }

		  splice_state = 0 ; //exon again
      }

    }
  
    // Right side not aligned
    for (int dna_pos=dstart+1; dna_pos<dna_len; dna_pos++) {
		s_align[dna_pos] = 4; 
    }
    for (int est_pos=rstart+1; est_pos<est_len; est_pos++) {
		e_align[est_pos] = 4;
    } 

	// for (int n=0;n<dna_len;n++)
	// 	fprintf(stdout,"%i",s_align[n]);
	// fprintf(stdout,"\n");
	// for (int n=0;n<est_len;n++)
	// 	fprintf(stdout,"%i",e_align[n]);
	// fprintf(stdout,"\n");
      
    return 0;
  }

  // No alignment found
  else{
	  for (int dna_pos=0; dna_pos<dna_len; dna_pos++) {
		  s_align[dna_pos] = 4 ; 
	  }
	  for (int est_pos=0; est_pos<est_len; est_pos++) {
		  e_align[est_pos] = 4 ;
	  }
	  //fprintf(stdout,"No alignment for this read\n");
	  return 1;
  }

}
