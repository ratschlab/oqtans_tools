#include <sys/time.h>
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include <string>
#include <time.h>
#include <stdlib.h>
#include <wait.h>

#include "palmapper.h"
#include "print.h"
#include "dyn_prog/qpalma_dp.h"

#include <palmapper/Genome.h>
#include <palmapper/Hits.h>
#include <palmapper/TopAlignments.h>
#include <palmapper/Util.h>



TopAlignments::TopAlignments(GenomeMaps *genomemaps_) :
	top_alignments(),
	num_spliced_alignments(0),
	num_unspliced_alignments(0),
	verbosity(0),
	num_filtered(0),
	current_ind(0),
	temp_ind(0)
{
	//int ret = pthread_mutex_init(&top_mutex, NULL) ;// = PTHREAD_MUTEX_INITIALIZER ;
	//assert(ret==0) ;
	
	ALIGNSEQ = (char *) malloc((_config.MAX_READ_LENGTH + 3 * Config::MAX_EDIT_OPS)
							   * sizeof(char));
	if (ALIGNSEQ == NULL) {
		fprintf(stderr, "[init_alignment_structures] Could not allocate memory\n");
		exit(1);
	}
	
	genomemaps = genomemaps_ ;
	
}

u_int8_t TopAlignments::report_unspliced_hit(Read const &read, HIT *hit, int num, QPalma const *qpalma)
{
	alignment_t *algn_hit = gen_alignment_from_hit(read, hit, qpalma) ;
	algn_hit->from_gm = 2;
	algn_hit->hit = hit ;
	algn_hit->num = num ;
	
	//printhit(read,hit);
	
	//if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP)
	
	add_alignment_record(algn_hit, 1) ;
	return 1 ;
}

void TopAlignments::determine_transcription_direction(char strand,char orientation, int side, char &transcription, char &read_forward)
{
	
	if (_config.PROTOCOL== PROTOCOL_UNSTRANDED){
		transcription=strand;
	}
	
	else{
		
		//Left reads
		if (side==1){
			//For sense transcription, the left read is in the same direction
			if (strand==orientation)
				transcription=(_config.PROTOCOL== PROTOCOL_FIRST)?'+':'-';
			else
				transcription=(_config.PROTOCOL==PROTOCOL_FIRST)?'-':'+';
		}
		else{
			//Right reads
			if (side==0){
				//For antisense transcription, the right read is in the same direction
				if (strand==orientation)
					transcription=(_config.PROTOCOL==PROTOCOL_FIRST)?'-':'+';
				else
					transcription=(_config.PROTOCOL==PROTOCOL_FIRST)?'+':'-';
			}
			//No information about the reads: the strand of a spliced alignment gives the direction of the transcription
			else
				transcription=strand;
		}
	}
	
	//For all cases, the reverse bit has to be given according to the coding strand: reverse the read orientation if found on the non coding strand
	if(strand=='-')
		read_forward=(orientation=='+')?'-':'+';
	else
		read_forward=orientation;

}



int TopAlignments::construct_aligned_string(Read const &read, HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *qual_mismatches_p, int *num_matches_p)
{
	int alignment_length = read.length() ;
	
	int j;
	int count_char = 0;
	char gap_offset = 0;
	char gap_in_read = 0;
	char gap_in_chr = 0;

	int num_gaps = 0 ;
	int num_mismatches = 0 ;
	int qual_mismatches = 0 ;
	int num_matches = read.length() ;

	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;

	if (hit->orientation == '+') {
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else {
		readstart = hit->start - (((int)read.length()) - hit->readpos - hitlength + 2)
			+ hit->start_offset; // 0-initialized
	}


	// #########################################
	// ###### PERFECT ALIGNMENTS: #############

	

	if (hit->mismatches == 0) {

		count_char = 0;
		
		for (size_t i=0; i < read.length(); i++) {
			char readbase;
			readbase = (hit->orientation == '+'? read.data()[i] : get_compl_base(read.data()[read.length()-i-1]));
			char readqual = (hit->orientation == '+'? read.quality(0)[i] : read.quality(0)[read.length()-i-1]);
			if ((*hit->chromosome)[readstart+i] != readbase) 
			{
				if (_config.BSSEQ) 
					sprintf(ALIGNSEQ + count_char, "{%c%c}", (*hit->chromosome)[readstart+i], readbase);
				else
					sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart+i], readbase);
				count_char += 4;
				if (readbase!='N' && (*hit->chromosome)[readstart+i]!='N')
				{
					num_mismatches += 1 ;
					qual_mismatches += readqual - read.get_quality_offset() ;
				}
				num_matches-- ;
			}
			else {
				ALIGNSEQ[count_char] = (*hit->chromosome)[readstart+i];
				count_char++;
			}
		}
		ALIGNSEQ[count_char] = '\0';
		


		if (num_mismatches_p)
			*num_mismatches_p = num_mismatches ;
		if (qual_mismatches_p)
			*qual_mismatches_p = qual_mismatches ;
		if (num_matches_p)
			*num_matches_p = num_matches ;
		
		
		return alignment_length ;
		
	}
	// #########################################
	// #########################################



	
	// fprintf(stdout, "CHR SEQ: ");
	// for (int i=0; i<(int)read.length()+5;i++)
	// 	fprintf(stdout, "%c",get_compl_base(hit->chromosome->operator []((int)read.length()+4+readstart-i)));
	// fprintf(stdout, "\n");




	// #########################################
	// ###### ALIGNMENTS WITH MISMATCHES: ######


	for (int i=0; i<hit->mismatches; i++)
		assert(hit->edit_op[i].pos>=-((int) read.length()) && hit->edit_op[i].pos<=((int)read.length())) ;
	
	// sort mismatches in ascending order according to their abs positions and 'gap before mm'-strategy if equal
	qsort(hit->edit_op, hit->mismatches, sizeof(EDIT_OPS), compare_editops);

	for (int i=0; i<hit->mismatches; i++)
		assert(hit->edit_op[i].pos>=-((int) read.length()) && hit->edit_op[i].pos<=((int)read.length())) ;

	ALIGNSEQ[0] = '\0';
	
	for (j = 0; j != hit->mismatches; ++j) 
    {
		//fprintf(stderr, "j=%i, edit_op[j].pos=%i\n count_char=%i", j, hit->edit_op[j].pos, count_char) ;
		assert(hit->edit_op[j].pos>= -((int)read.length()) && hit->edit_op[j].pos<=((int)read.length())) ;
		
		if (hit->edit_op[j].pos < 0) 
		{
			hit->edit_op[j].pos = -hit->edit_op[j].pos;
			gap_in_chr = 1;
		}

		if (j == 0) 
		{
			if (hit->edit_op[0].pos - 1>=0)
			{	// from start to first mismatch:
				for (size_t i=0; (int)i<hit->edit_op[0].pos - 1; i++) {
					char readbase;
					readbase = (hit->orientation == '+'? read.data()[i] : get_compl_base(read.data()[read.length()-i-1]));
					char readqual = (hit->orientation == '+'? read.quality(0)[i] : read.quality(0)[read.length()-i-1]);
					if ((*hit->chromosome)[readstart+i] != readbase) {
						if (_config.BSSEQ)
							sprintf(ALIGNSEQ + count_char, "{%c%c}", (*hit->chromosome)[readstart+i], readbase);
						else
							sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart+i], readbase);
						num_matches--;
						count_char+=4;
						if (readbase!='N' && (*hit->chromosome)[readstart +i]!='N'){
							num_mismatches ++;
							qual_mismatches += readqual - read.get_quality_offset() ;
						}
						
					}
					else{
						ALIGNSEQ[count_char] = (*hit->chromosome)[readstart+i] ;
						count_char++;						
					}
				}		   

			}
		} 
		else if (hit->edit_op[j].pos - hit->edit_op[j - 1].pos != 0) {	// from mismatch to mismatch
			for (size_t i=0; (int)i<hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read; i++) {
				char readbase;
				readbase = (hit->orientation == '+'? read.data()[i + hit->edit_op[j-1].pos - gap_in_read] :
							get_compl_base(read.data()[read.length() - hit->edit_op[j-1].pos - i - 1 + gap_in_read]));
				char readqual = (hit->orientation == '+'? read.quality(0)[i + hit->edit_op[j-1].pos - gap_in_read] : read.quality(0)[read.length() - hit->edit_op[j-1].pos - i - 1 + gap_in_read]);
				if ((*hit->chromosome)[readstart+i+hit->edit_op[j-1].pos+gap_offset-gap_in_read] != readbase) {
					if (_config.BSSEQ)
						sprintf(ALIGNSEQ + count_char,"{%c%c}",(*hit->chromosome)[readstart+i+hit->edit_op[j-1].pos+gap_offset-gap_in_read],readbase);
					else
						sprintf(ALIGNSEQ + count_char,"[%c%c]",(*hit->chromosome)[readstart+i+hit->edit_op[j-1].pos+gap_offset-gap_in_read],readbase);
					num_matches--;					
					count_char += 4;
					if ( readbase!='N' && (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read)+i]!='N' )	{						
						num_mismatches++;	
						qual_mismatches += readqual - read.get_quality_offset() ;
					}
					
					
				}
				else {
					ALIGNSEQ[count_char] = readbase;
					count_char++;
				}
			}

		} // else: edit_op[j-1] must have been a gap!

		gap_in_read = 0;

		if (hit->edit_op[j].mm)
		{
			num_matches-- ;
			
			if (hit->orientation == '+')
			{
				sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset], read.data()[hit->edit_op[j].pos - 1]);
				if ( read.data()[hit->edit_op[j].pos - 1]!='N' && (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset]!='N' )
				{
					num_mismatches++ ;
					qual_mismatches += read.quality(0)[hit->edit_op[j].pos - 1] - read.get_quality_offset() ;
			   
				}
			}
			else
			{
			
				sprintf(ALIGNSEQ + count_char, "[%c%c]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset], get_compl_base(read.data()[((int)read.length()) - hit->edit_op[j].pos]));
				if ( (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset]!='N' && get_compl_base(read.data()[((int)read.length()) - hit->edit_op[j].pos ])!='N' ){					
					num_mismatches++ ;
					qual_mismatches += read.quality(0)[((int)read.length()) - hit->edit_op[j].pos] - read.get_quality_offset() ;
				}	
			}
			
		}
		else if (gap_in_chr)
		{
			num_gaps++ ;
			num_matches-- ;
			alignment_length-- ;

			if (hit->orientation == '+')
			{
				if (read.data()[hit->edit_op[j].pos - 1]!='N'){					
					qual_mismatches += read.quality(0)[hit->edit_op[j].pos - 1] - read.get_quality_offset() ;
				}
				
			}
			else
			{
				if (read.data()[read.length() - hit->edit_op[j].pos]!='N'){					
					qual_mismatches += read.quality(0)[read.length() - hit->edit_op[j].pos] - read.get_quality_offset() ;
				}
				
			}
			
			if (_config.BSSEQ && (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset] == (hit->orientation == '+'? read.data()[hit->edit_op[j].pos - 1] : get_compl_base(read.data()[read.length() - hit->edit_op[j].pos])))
			{
				// convert [-C]{CT} to C[-T] and [-G]{GA} to G[-A]:
				sprintf(ALIGNSEQ + count_char, "%c[-%c]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset], (hit->orientation == '+' ? read.data()[hit->edit_op[j].pos] : get_compl_base(read.data()[read.length() - hit->edit_op[j].pos - 1])));
				count_char++;
				hit->edit_op[j].pos++;
			}
			else
			{
				if (hit->orientation == '+')
				{
					sprintf(ALIGNSEQ + count_char, "[-%c]", read.data()[hit->edit_op[j].pos - 1]);
				}
				else
				{
					sprintf(ALIGNSEQ + count_char, "[-%c]", get_compl_base(read.data()[read.length() - hit->edit_op[j].pos]));
				}
			}
			
			gap_offset--;
			gap_in_chr = 0;
		} else
		{
			num_gaps++ ;
			alignment_length++ ;
			
			sprintf(ALIGNSEQ + count_char, "[%c-]", (*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset]);

			gap_offset++;
			gap_in_read = 1;
		}

		count_char += 4;
    }

	// from last mismatch to end of read:

	if (((int)read.length()) - hit->edit_op[j - 1].pos + gap_in_read >= 0){
		for (size_t i=0; (int)i < ((int) read.length())  - hit->edit_op[j - 1].pos + gap_in_read; i++) {
			char readbase;
			readbase = (hit->orientation == '+'?
						read.data()[i + hit->edit_op[j-1].pos - gap_in_read] :
						get_compl_base(read.data()[read.length() - hit->edit_op[j-1].pos - i - 1 + gap_in_read]));
			char readqual = (hit->orientation == '+'? read.quality(0)[i + hit->edit_op[j-1].pos - gap_in_read] : read.quality(0)[read.length() - hit->edit_op[j-1].pos - i - 1 + gap_in_read]);
			if ((*hit->chromosome)[readstart+i+hit->edit_op[j-1].pos+gap_offset-gap_in_read] != readbase) {
				if (_config.BSSEQ) 
					sprintf(ALIGNSEQ + count_char,"{%c%c}",	(*hit->chromosome)[readstart+i+hit->edit_op[j-1].pos+gap_offset-gap_in_read],readbase);
				else
					sprintf(ALIGNSEQ + count_char,"[%c%c]",	(*hit->chromosome)[readstart+i+hit->edit_op[j-1].pos+gap_offset-gap_in_read],readbase);
				count_char += 4;
				num_matches--;
				if ( readbase!='N' && (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read)+i]!='N' ){
					num_mismatches++;
					qual_mismatches += readqual - read.get_quality_offset() ;					
				}
				
			}
			else {
				ALIGNSEQ[count_char] = readbase;
				count_char++;
			}
		}
	}
	
	
	ALIGNSEQ[count_char] = '\0';
//	fprintf(stdout, "start=%i hit_len=%i count_char=%i hit_mm=%i num_mismatches=%i num_matches=%i\n", hit->start, hitlength, count_char, hit->mismatches, num_mismatches, num_matches) ;
//	fprintf(stdout,"%s\n",ALIGNSEQ);

	// #########################################
	// #########################################

	
	if (num_gaps_p)
		*num_gaps_p = num_gaps ;
	if (num_mismatches_p)
		*num_mismatches_p = num_mismatches ;
	if (qual_mismatches_p)
		*qual_mismatches_p = qual_mismatches ;
	if (num_matches_p)
		*num_matches_p = num_matches ;
	

	
	return alignment_length ;

}

alignment_t *TopAlignments::gen_alignment_from_hit(Read const &read, HIT *best_hit, QPalma const * qpalma)
{
	int hitlength = best_hit->end - best_hit->start + 1;
	unsigned int readstart;
	unsigned int readend;

	int num_gaps=0, num_mismatches=0, num_matches=((int)read.length()), qual_mismatches=0 ;
	int alignment_length=construct_aligned_string(read, best_hit, &num_gaps, &num_mismatches, &qual_mismatches, &num_matches); // constructs the aligned string with mismatches and gaps in ALIGNSEQ

	if (best_hit->orientation == '+') 
    {
		readstart = best_hit->start - best_hit->readpos + best_hit->start_offset; // start pos of read in genome	0-initialized
		readend = readstart + alignment_length ;
    } else 
    {
		readstart = best_hit->start - (((int)read.length()) - best_hit->readpos - hitlength + 2) + best_hit->start_offset; // 0-initialized
		readend = readstart + alignment_length ;
    }

	
	alignment_t *best = new alignment_t ;
	best->qpalma_score = 1000 ;
	best->sort_key = random()/RAND_MAX ;
	best->num_matches = num_matches ;
	best->num_mismatches = num_mismatches ;
	best->qual_mismatches = qual_mismatches ;
	best->num_gaps = num_gaps ;
	best->chromosome = best_hit->chromosome ;
	best->orientation = best_hit->orientation ;
	best->strand='+' ;
	strcpy(best->read_id, read.id()) ;
	best->exons.push_back(readstart) ;
	best->exons.push_back(readend) ;
	strcpy(best->read_anno, ALIGNSEQ) ;
	best->min_exon_len=best->exons[1]-best->exons[0] ;
	

	if (_config._personality == Palmapper) {
		best->max_intron_len = 0 ;
		best->spliced = false ;
		best->non_consensus_alignment = false ;
		best->remapped=false ;

		bool alignment_passed_filters= (num_mismatches <= _config.NUM_MISMATCHES && num_gaps <= _config.NUM_GAPS && num_mismatches+num_gaps <= _config.NUM_EDIT_OPS) ;
		best->passed_filters=alignment_passed_filters ;

		if (qpalma){
			double score1=qpalma->score_unspliced(read, ALIGNSEQ,'+',best->orientation) ;
			double score2=qpalma->score_unspliced(read, ALIGNSEQ,'-',best->orientation) ;
			
			if (score1>score2)
				best->qpalma_score =score1;
			else{
				best->qpalma_score =score2;
				best->strand='-';
				best->orientation=(best->orientation=='+')?'-':'+';
			}
			
		}
		
		else
			best->qpalma_score = - ((best_hit->mismatches * _config.MM_SCORE) + (best_hit->gaps * _config.GAP_SCORE) + ((hitlength-best_hit->mismatches-best_hit->gaps) * _config.M_SCORE));
	} else {
		best->qpalma_score = - ((best_hit->mismatches * _config.MM_SCORE) + (best_hit->gaps * _config.GAP_SCORE) + ((hitlength-best_hit->mismatches-best_hit->gaps) * _config.M_SCORE));
	}
	best->rtrim_cut=0 ;
	best->polytrim_cut_start=0 ;
	best->polytrim_cut_end=0 ;

	best->hit = NULL ;
	best->num = 0 ;

	return best ;
}

// returns > 0 if a1 scores better than a2
// returns < 0 if a1 scores worse than a2
// returns 0 if a1 and a2 score equally well
int32_t TopAlignments::compare_score(alignment_t *a1, alignment_t *a2) {

	assert(a1->qpalma_score!=1000) ;
	assert(a2->qpalma_score!=1000) ;

	if (fabs(a1->qpalma_score-a2->qpalma_score)<1e-6)

//#else //TODO joerg is it really important to use integer abs here when using personality genomemapper?
//	if (abs(a1->qpalma_score-a2->qpalma_score)<1e-6)
//#endif
    {
		// scores are numerically identical, use start position as sorting criterion
		if (a1->exons[0] > a2->exons[0])
			return 1;
		else if (a1->exons[0] < a2->exons[0])
			return -1;
		else
			return 0;
    }
	
	
	if (a1->qpalma_score > a2->qpalma_score)
		return 1;
	else if (a1->qpalma_score < a2->qpalma_score)
		return -1;
	else
		// a1->score == a2->score
		return 0;
}



int TopAlignments::spliced_is_overlapping(alignment_t *a1, alignment_t *a2) 
{

	if (!a1->spliced || !a2->spliced)
		return 0;
	if(a1->non_consensus_alignment!=a2->non_consensus_alignment)
		return 0;
	if(a1->chromosome!=a2->chromosome)
		return 0;
	
	int startalign1=a1->exons[0];
	int endalign1=a1->exons[a1->exons.size()-1];
	int startalign2=a2->exons[0];
	int endalign2=a2->exons[a2->exons.size()-1];
	if (endalign1<startalign2 || startalign1 > endalign2)	
		return 0;
	
	if (a2->passed_filters && !a1->passed_filters)
		return 1;

	if (a1->passed_filters && !a2->passed_filters)
		return -1;

	if(a1->qpalma_score < a2->qpalma_score) 
		return 1;
	
	return -1;
}


int TopAlignments::alignment_is_opposite(alignment_t *a1, alignment_t *a2) 
{
	if (a1->num_matches!=a2->num_matches) 
		return 0 ;

	if (a1->strand==a2->strand && a1->orientation==a2->orientation)
		return 0 ;
	if (a1->chromosome!=a2->chromosome)
		return 0 ;
	if (a1->exons!=a2->exons)
		return 0 ;

	if (a1->non_consensus_alignment && !a2->non_consensus_alignment)
		return 1 ;
	if (!a1->non_consensus_alignment && a2->non_consensus_alignment)
		return -1 ;

	if(a1->qpalma_score < a2->qpalma_score) 
		return 1;

	//if (a1->strand=='+') // ???
	//	return -1 ;
	
	return -1 ;
}

bool TopAlignments::alignment_is_equal(alignment_t *a1, alignment_t *a2) 
{
	
	
	/*if (fabs(a1->qpalma_score-a2->qpalma_score)>1e-3){
		return false ;
		}*/ // cant use this due to the remapping scoring
	
	if (a1->num_matches!=a2->num_matches){
		return false ;
	}
	if (a1->orientation!=a2->orientation){
		return false ;
	}
	if (a1->strand!=a2->strand){
		return false ;
	}
	if (a1->chromosome!=a2->chromosome){
		return false ;
	}
	if (a1->exons!=a2->exons)
    {
		// for (int i=0; i<a1->exons.size(); i++)
		// 	fprintf(stderr, "a1.e[%i]=%i\n", i, a1->exons[i]) ;
		// for (int i=0; i<a2->exons.size(); i++)
		// 	fprintf(stderr, "a2.e[%i]=%i\n", i, a2->exons[i]) ;
		
		return false ;
    }

	return true ;	
}

int TopAlignments::non_consensus_overlaps_consensus(alignment_t *a1, alignment_t *a2)
{


	if (a1->non_consensus_alignment == a2->non_consensus_alignment)
		return 0;
	//if (a1->strand==a2->strand || a1->orientation != a2->orientation)
	//	return 0;
	if (a1->chromosome != a2->chromosome)
		return 0;
	if (!a1->spliced && !a2->spliced)
		return 0;
	
	int startalign1=a1->exons[0];
	int endalign1=a1->exons[a1->exons.size()-1];
	int startalign2=a2->exons[0];
	int endalign2=a2->exons[a2->exons.size()-1];

	//fprintf(stdout,"a1: %i-%i\n",startalign1,endalign1);
	//fprintf(stdout,"a2: %i-%i\n",startalign2,endalign2);

	if (endalign1<startalign2 || startalign1 > endalign2)
		return 0;
	
	if (a1->non_consensus_alignment)
		return 1;
	
	return -1;
}


void TopAlignments::init_top_alignment_indice() 
{
	current_ind=0;
	temp_ind=0;
	
}

void TopAlignments::update_top_alignment_indice() 
{
	//fprintf(stdout,"temp indice %i\n",temp_ind);
	
	current_ind=temp_ind;
	
}

void TopAlignments::clean_top_alignment_record() 
{
	//pthread_mutex_lock( &top_mutex) ;
	
	// Cleaning up and starting over.
	for (uint32_t i = 0; i < top_alignments.size(); i++) 
    {
		delete top_alignments[i];
    }

	top_alignments.clear();
	num_spliced_alignments = 0;
	num_unspliced_alignments = 0;

	num_filtered=0;
	current_ind=0;
	temp_ind=0;

	//pthread_mutex_unlock( &top_mutex) ;
}

void TopAlignments::start_top_alignment_record() 
{
	clean_top_alignment_record() ;
}

void TopAlignments::check_alignment(struct alignment_t * alignment)
{
	//if (alignment->exons.size()<4)
	//	return ;
	
	for (int i=0; i<((int)alignment->exons.size())-1; i++)
    {
		assert(alignment->exons[i]<alignment->exons[i+1]) ;
    }

}


/** performs a quicksort on an array output of length size
 * it is sorted from in ascending (for type T) */
void TopAlignments::qsort_top_alignments(alignment_t** output, int num) {
	
	if(num == 1)
		return;
	

	if (num == 2) {
		if (compare_score(output[0],output[1])<0) {
			alignment_t *c = output[0];
			output[0] = output[1];
			output[1] = c;
		}
		return;
	}
	alignment_t *split = output[num / 2];

	int32_t left = 0;
	int32_t right = num - 1;

	while (left <= right) {
		while (compare_score(output[left],split)>0)
			left++;
		while (compare_score(output[right],split)<0)
			right--;

		if (left <= right) {
			alignment_t *c = output[left];
			output[left] = output[right];
			output[right] = c;
			left++;
			right--;
		}
	}

	if (right + 1 > 1)
		qsort_top_alignments(output, right + 1);

	if (num - left > 1)
		qsort_top_alignments(&output[left], num - left);
}

void TopAlignments::sort_top_alignment_list()
{
	
	//fprintf(stdout, "Number of alignments %i\n", top_alignments.size());

	// Remove unfiltered alignments from non_consensus_search
	if (_config.non_consensus_search)
	{
		for (size_t i=0; i<top_alignments.size();++i){
			if (!top_alignments[i]->passed_filters){
				delete top_alignments[i] ;
				top_alignments[i]=NULL ;
				top_alignments.erase(top_alignments.begin()+i) ;
				i--;
			}
		}
		
	}

	// Sort by qpalma score and positions remaining alignments if not all reported
	{
		
		alignment_t ** arr = NULL ;
		size_t nbr_aligments=top_alignments.size();
		try 
		{
			arr = new alignment_t*[nbr_aligments] ;			
		}
		catch (std::bad_alloc&)
		{
			fprintf(stderr, "[sort_top_alignment_list] ERROR Could not allocate memory\n");
		}
		
		for (int i = 0; i < (int)nbr_aligments; i++)
		{
			arr[i] = top_alignments[i];
			if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP || _config.OUTPUT_FILTER==OUTPUT_FILTER_ALL)
				arr[i]->sort_key = arr[i]->qpalma_score ;
			if (_config.OUTPUT_FILTER==OUTPUT_FILTER_LIMIT)
				arr[i]->sort_key = random()/RAND_MAX ;
		}
		
		qsort_top_alignments(arr, nbr_aligments);
		
		for (int i = 0; i < (int)nbr_aligments; i++)
			top_alignments[i] = arr[i];
		
		delete[] arr;
	}
	
	if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP || _config.OUTPUT_FILTER==OUTPUT_FILTER_LIMIT)
	{
		unsigned int limit = _config.OUTPUT_FILTER_NUM_TOP ;
		
		if (_config.OUTPUT_FILTER==OUTPUT_FILTER_LIMIT)
			limit = _config.OUTPUT_FILTER_NUM_LIMIT ;
		
		if (top_alignments.size() > limit)
		{
			for(size_t i=top_alignments.size()-1; i>=limit; i--)
			{
				delete top_alignments[i] ;
				top_alignments[i]=NULL ;
				top_alignments.pop_back();
			}
		}
		assert(top_alignments.size() <= limit) ;
	}
}

void TopAlignments::end_top_alignment_record(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP, int rtrim_cut, int polytrim_cut_start, int polytrim_cut_end, JunctionMap &junctionmap) {

	if (top_alignments.empty() && ! _config.INCLUDE_UNMAPPED_READS_SAM)
		return;

	//pthread_mutex_lock( &top_mutex) ;

	sort_top_alignment_list();
	
	for (unsigned int i=0; i<top_alignments.size(); i++)
    {
		top_alignments[i]->rtrim_cut=rtrim_cut ;
		top_alignments[i]->polytrim_cut_start=polytrim_cut_start ;
		top_alignments[i]->polytrim_cut_end=polytrim_cut_end ;
		if (_config._personality == Palmapper)
			check_alignment(top_alignments[i]) ;
	}

	if (_config._personality == Palmapper)
	{
		if (_config.REPORT_SPLICED_READS)
		{
			for (unsigned int i=0; i<top_alignments.size(); i++)
				if (top_alignments[i]->exons.size()>2)
					genomemaps->report_spliced_read(*top_alignments[i]->chromosome, top_alignments[i]->exons, 
													top_alignments[i]->num_matches, i) ;
		}
		if (_config.REPORT_MAPPED_READS)
		{
			for (unsigned int i=0; i<top_alignments.size(); i++)
			{
				if (top_alignments[i]->exons.size()<=2) 
					genomemaps->report_mapped_read(*top_alignments[i]->chromosome, top_alignments[i]->exons[0], top_alignments[i]->exons[1], 
												   top_alignments[i]->num_matches, i) ;
			}
		}
		
		if (_config.REPORT_JUNCTIONS)
		{
			for (unsigned int i=0; i < top_alignments.size() && i < 1; i++)
			{
				if (top_alignments[i]->remapped) // only count real alignments, not remapped alignments
					continue ;
				for (unsigned int j=2; j < top_alignments[i]->exons.size(); j+=2 ){
					junctionmap.insert_junction(top_alignments[i]->strand,top_alignments[i]->chromosome->nr(), 
												top_alignments[i]->exons[j-1], top_alignments[i]->exons[j]-1,
												!top_alignments[i]->non_consensus_intron[j/2-1], 
												top_alignments[i]->intron_consensus[j/2-1].c_str(), 
												top_alignments[i]->min_exon_len, 
												top_alignments[i]->read_id, 1) ;
				}
			}
		}
    }

	print_top_alignment_records(read, OUT_FP, SP_OUT_FP) ;
	
	//pthread_mutex_unlock( &top_mutex) ;

	clean_top_alignment_record();
}

alignment_t * TopAlignments::add_alignment_record(alignment_t *alignment, int num_alignments)
{

	std::vector<alignment_t *>::iterator it;

	if (alignment == NULL)
		return NULL;
	assert(num_alignments>0);

	if (_config._personality != Palmapper || !alignment->spliced)
		num_unspliced_alignments += num_alignments;
	else
		num_spliced_alignments += num_alignments;


	if (_config._personality == Palmapper) {
		// if (!alignment->passed_filters)
		// {
		// 	delete alignment;
		// 	return NULL ;
		// }
		check_alignment(alignment) ;
	}


	if (verbosity >= 2)
		printf("entering alignment with score %f\n", alignment->qpalma_score);
 
	// Go through the list and see whether we can find a hit that has a worse
	// score than this one.

	// Special case: list of top hits is still empty. Kick-start it with the
	// current hit.

	// fprintf(stdout,"number of alignments %i\n", (int)top_alignments.size());
	
	// fprintf(stdout,"#alignment (%i) with %i exons found for %s (score=%1.3f  matches=%i mismatches=%i gaps=%i strand=%c orientation=%c spliced=%i): %s\n",
	// 		alignment->non_consensus,
	// 		(int) alignment->exons.size() / 2, 
	// 		alignment->read_id, 
	// 		alignment->qpalma_score,
	// 		alignment->num_matches, 
	// 		alignment->num_mismatches, 
	// 		alignment->num_gaps, 
	// 		alignment->strand, 
	// 		alignment->orientation, 
	// 		alignment->spliced, 
	// 		alignment->read_anno);
	// for (size_t j = 0; j < alignment->exons.size(); j += 2)
	// 	fprintf(stdout, "# exon %i: %i - %i (%i)\n", (int)j / 2, alignment->exons[j], alignment->exons[j+ 1],alignment->chromosome->nr());

	// First alignment added
	if (top_alignments.empty())
    {
		top_alignments.push_back(alignment);
		current_ind=0;
		if (alignment->passed_filters)
			num_filtered++;
		return alignment ;
    }
	else
    {
		int chr = alignment->chromosome->nr();
		int endalign=alignment->exons[alignment->exons.size()-1];

		int current_chr=top_alignments[current_ind]->chromosome->nr();
		int current_start=top_alignments[current_ind]->exons[0];
		

		it = top_alignments.begin();

		// Sort alignment by chromosome and start position

		//Nothing to compare: add before
		if ( (chr < current_chr) || (chr==current_chr && endalign < current_start))
		{
			//fprintf(stdout,"   Add first\n");
			top_alignments.insert(it+current_ind,1,alignment);
			if (alignment->passed_filters)
				num_filtered++;
			return alignment;
		}
		else{
			
			int i =current_ind;
			//fprintf(stdout,"   Start comparing with %i\n", i);
			while (i< (int)top_alignments.size())
			{
				//fprintf(stdout,"   Compare with %i\n", i);
				// for (size_t j = 0; j < top_alignments[i]->exons.size(); j += 2)
				// 	fprintf(stdout, "# exon %i: %i - %i (%i)\n", (int)j / 2,  top_alignments[i]->exons[j],  top_alignments[i]->exons[j+ 1], top_alignments[i]->chromosome->nr());
				//Adapt pointer
				if (chr!=current_chr && chr==(int)top_alignments[i]->chromosome->nr()){
					current_chr = chr;
					current_ind = i;
				}
				
		
				//Alignment to compare is further away: no conflict possible
				if ((chr==(int)top_alignments[i]->chromosome->nr() && endalign < top_alignments[i]->exons[0]) 
					|| (int)top_alignments[i]->chromosome->nr() > chr)
				{
					//fprintf(stdout,"   Add before because no conflict possible\n");
					top_alignments.insert(it+i,1,alignment);
					if (temp_ind<i)
						temp_ind=i;
					else
						temp_ind++;
					
					if (alignment->passed_filters)
						num_filtered++;
					return alignment;
				}
	
					
					
				int ret;
				// Alignment already present
				if (alignment_is_equal(alignment,top_alignments[i]))
				{
					//fprintf(stdout,"   Ignore alignment: already present\n");
					if (temp_ind<i)
						temp_ind=i;
					
					if (alignment!=top_alignments[i])
						delete alignment ;
					return NULL ;
				}

				if (false)
				{
					//Overlapping consensus alignments: keep the one which passes filters with best score
					ret=spliced_is_overlapping(alignment, top_alignments[i]);
					if (ret==-1)
					{
						//fprintf(stdout,"   Replace alignment: spliced overlapping\n");
						if (!top_alignments[i]->passed_filters && alignment->passed_filters)
							num_filtered++;
						delete top_alignments[i] ;
						top_alignments[i]=alignment;
						if (temp_ind<i)
							temp_ind=i;
						
						return alignment;		
					}
					if (ret==1)
					{
						//fprintf(stdout,"   Ignore alignment: spliced overlapping\n");
						if (temp_ind<i)
							temp_ind=i;
						
						delete alignment;
						return NULL;						
					}
				}
				
				//Overlapping consensus alignment with a non consensus one: keep consensus
				ret=non_consensus_overlaps_consensus(alignment, top_alignments[i]);
				
				if (ret==-1)
				{
					//fprintf(stdout,"   Replace alignment: consensus and non consensus\n");
					if (!top_alignments[i]->passed_filters && alignment->passed_filters)
						num_filtered++;
					delete top_alignments[i] ;
					top_alignments[i]=alignment;
					if (temp_ind<i)
						temp_ind=i;
					
					return alignment;		
				}
				if (ret==1){
					//fprintf(stdout,"   Ignore alignment: consensus and non consensus\n");
					if (temp_ind<i)
						temp_ind=i;
					
					delete alignment;
					return NULL;						
				}
				
				// Opposite alignments:
				ret=alignment_is_opposite(alignment,top_alignments[i]);
				if (ret==-1)
				{
					//fprintf(stdout,"   Replace alignment: opposite\n");
					if (!top_alignments[i]->passed_filters && alignment->passed_filters)
						num_filtered++;
					delete top_alignments[i] ;
					top_alignments[i]=alignment;
					if (temp_ind<i)
						temp_ind=i;

					
					return alignment;		
				}
				if (ret==1)
				{
					//fprintf(stdout,"   Ignore alignment: opposite\n");
					if (temp_ind<i)
						temp_ind=i;
					
					delete alignment;
					return NULL;						
				}
			
				
				i++;
					
			}//End while i
				
			//Alignment is without conflict and from a further chromosome or position compared to the others
			//fprintf(stdout,"   Add at the end\n");
			top_alignments.push_back(alignment);
			temp_ind=top_alignments.size()-1;
			
			if (alignment->passed_filters)
				num_filtered++;
			return alignment;

		}
	}// Top alignment not empty

}


		   
int TopAlignments::print_top_alignment_records(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP)
{
	if (_config._personality == Palmapper && _config.OUTPUT_FORMAT==OUTPUT_FORMAT_BEDX)
    {
		return print_top_alignment_records_bedx(read, OUT_FP, SP_OUT_FP) ;
    }
	if (_config.OUTPUT_FORMAT==OUTPUT_FORMAT_SHORE || _config.OUTPUT_FORMAT==OUTPUT_FORMAT_BED)
    {
		return print_top_alignment_records_shorebed(read, OUT_FP, SP_OUT_FP) ;
    }

	if ( _config.OUTPUT_FORMAT==OUTPUT_FORMAT_SAM || _config.OUTPUT_FORMAT==OUTPUT_FORMAT_BAM )
    {
		return print_top_alignment_records_sam(read, OUT_FP, SP_OUT_FP) ;
    }

	
	fprintf(stderr, "ERROR: unknown output format\n") ;
	exit(1) ;
	return 0 ;
}

int TopAlignments::print_top_alignment_records_bedx(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP)
{

	int top_unspliced_alignments=0;
	int top_spliced_alignments=0;

	if (top_alignments.size()==0)
		return 0 ;
	alignment_t * best = top_alignments[0] ;
	assert(best->exons.size() >= 2);

	if (best->rtrim_cut!=0)
		assert(best->polytrim_cut_start==0 && best->polytrim_cut_end==0) ;
	if (best->polytrim_cut_start!=0)
		assert(best->polytrim_cut_end==0) ;

	std::ostream* MY_OUT_FP = OUT_FP ;

	pthread_mutex_lock( &_stats.alignment_num_mutex ) ;

	if (best->spliced)
    {
		assert(_config.SPLICED_HITS) ;
		MY_OUT_FP = SP_OUT_FP ;
		_stats.alignment_num_spliced_best++ ;
		++top_spliced_alignments;
		
    }
	else{
		
		_stats.alignment_num_unspliced_best++ ;
		++top_unspliced_alignments;
	}
	

	for (int i=1; i<(int)top_alignments.size(); i++)
		if (top_alignments[i]->spliced){
			++top_spliced_alignments;
			_stats.alignment_num_spliced_suboptimal+=1 ;
		}
		else{
			_stats.alignment_num_unspliced_suboptimal+= 1 ;
			++top_unspliced_alignments;
		}
	
	_stats.print_alignment_stats() ;
	
	pthread_mutex_unlock( &_stats.alignment_num_mutex) ;

	char transcription_direction;
	char read_orientation;
	determine_transcription_direction(best->strand,best->orientation,_config.STRAND,transcription_direction,read_orientation);

	// Print spliced alignment to open file MY_OUT_FP in BED format
	fprintf(MY_OUT_FP, "%s\t%s\t%d\t%d\t%s\t%d\t%c\t%i\t%i\t0,0,0\t%d\t",
			best->from_gm == 2? "gm":"qp",
			best->chromosome->desc(),
			best->exons[0],
			best->exons[best->exons.size()-1],
			best->read_id,
			best->num_matches,
			transcription_direction,
			top_unspliced_alignments,
			top_spliced_alignments,
			int(best->exons.size() / 2));



	
	fprintf(MY_OUT_FP, "%d", best->exons[1] - best->exons[0]);
	assert(best->exons[1] - best->exons[0] >  0) ;

	
	for (uint32_t i = 2; i < best->exons.size(); i+=2)
    {
		fprintf(MY_OUT_FP, ",%d", best->exons[i+1] - best->exons[i]);
		assert(best->exons[i+1] - best->exons[i] >  0) ;

    }

	fprintf(MY_OUT_FP, "\t0");
	for (uint32_t i = 2; i < best->exons.size(); i+=2)
    {
		fprintf(MY_OUT_FP, ",%d", best->exons[i] - best->exons[0]);
		assert(best->exons[i] - best->exons[0] > 0) ;
    }
	
	{
		int polytrim_cut_start=best->polytrim_cut_start ;
		int polytrim_cut_end=best->polytrim_cut_end ;

		if (_config.RTRIM_STRATEGY && !_config.POLYTRIM_STRATEGY)
		{
			assert(polytrim_cut_end==0) ;
			polytrim_cut_end = best->rtrim_cut ;
		}

		char *read_anno=new char[strlen(best->read_anno)+(polytrim_cut_start+polytrim_cut_end)*4+20] ;
		const char *read_qual = read.quality(0);
		int read_len = read.length() ;
		
		if ((_config.POLYTRIM_STRATEGY||_config.RTRIM_STRATEGY) && (polytrim_cut_start>0 || polytrim_cut_end>0))
		{
			char const * orig_read = read.get_orig()->data() ;
			int orig_len = read.get_orig()->length() ;
			read_len = orig_len ;
			read_qual = read.get_orig()->quality(0);

			if (read_orientation=='+')
			{
				for (int i=0; i<polytrim_cut_start; i++)
				{
					read_anno[4*i]='[' ;
					read_anno[4*i+1]='=' ;
					read_anno[4*i+2]=orig_read[i] ;
					read_anno[4*i+3]=']' ;
				}
				strcpy(&(read_anno[4*polytrim_cut_start]), best->read_anno) ;
				int len=strlen(read_anno) ;
				for (int i=0; i<polytrim_cut_end; i++)
				{
					read_anno[len+4*i]='[' ;
					read_anno[len+4*i+1]='=' ;
					read_anno[len+4*i+2]=orig_read[orig_len-i-1] ;
					read_anno[len+4*i+3]=']' ;
				}
				read_anno[len+4*polytrim_cut_end]=0 ;
			}
			else
			{ // in this case read_anno is reverse complemented
				for (int i=0; i<polytrim_cut_end; i++)
				{
					read_anno[4*i]='[' ;
					read_anno[4*i+1]='=' ;
					read_anno[4*i+2]=get_compl_base(orig_read[orig_len-i-1]) ;
					read_anno[4*i+3]=']' ;
				}
				strcpy(&read_anno[4*polytrim_cut_end], best->read_anno) ;
				int len=strlen(read_anno) ;
				for (int i=0; i<polytrim_cut_start; i++)
				{
					read_anno[len+4*i]='[' ;
					read_anno[len+4*i+1]='=' ;
					read_anno[len+4*i+2]=get_compl_base(orig_read[i]) ;
					read_anno[len+4*i+3]=']' ;
				}
				read_anno[len+4*polytrim_cut_start]=0 ;
			}
		}
		else
			strcpy(read_anno, best->read_anno) ;
		
		double qpalma_score = best->qpalma_score ;
		
		{
			if (read_orientation=='+')
				fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s;quality=%s", 
						qpalma_score, best->num_matches, best->num_gaps, best->min_exon_len, best->max_intron_len, read_orientation, read_anno, read_qual) ;
			else
			{
				// reverse order of quality 
				char qual[read_len+1] ;
				for (int i=0; i<read_len; i++)
					qual[i]=read_qual[((int)read_len)-i-1] ;
				qual[read_len]=0 ;
				
				fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s;quality=%s", 
						qpalma_score, best->num_matches, best->num_gaps, best->min_exon_len, best->max_intron_len, read_orientation, read_anno, qual) ;
			}
			if (_config.POLYTRIM_STRATEGY)
			{
				if (best->polytrim_cut_start)
					fprintf(MY_OUT_FP, ";polytrimStart=%i", best->polytrim_cut_start) ;
				if (best->polytrim_cut_end)
					fprintf(MY_OUT_FP, ";polytrimEnd=%i", best->polytrim_cut_end) ;
			}
			if (_config.RTRIM_STRATEGY)
			{
				if (best->rtrim_cut)
					fprintf(MY_OUT_FP, ";rtrimEnd=%i", best->rtrim_cut) ;
			} 
		}
		delete[] read_anno ;
	}


	for (unsigned int j=1; j<top_alignments.size(); j++)
    {
		alignment_t * second  = top_alignments[j] ;
		
		assert(second->exons.size() >= 2);

		determine_transcription_direction(second->strand,second->orientation,_config.STRAND,transcription_direction,read_orientation);

		fprintf(MY_OUT_FP, "\t%s\t%s\t%d\t%d\t%d\t%c\t%d\t",
				second->from_gm == 2? "gm":"qp",
				second->chromosome->desc(),
				second->exons[0],
				second->exons[second->exons.size() - 1],
				second->num_matches,
				transcription_direction,
				int(second->exons.size() / 2));



		fprintf(MY_OUT_FP, "%d", second->exons[1] - second->exons[0]);
		assert(second->exons[1] - second->exons[0] >  0) ;

		
		for (uint32_t i = 2; i < second->exons.size(); i+=2)
		{
			fprintf(MY_OUT_FP, ",%d", second->exons[i+1] - second->exons[i]);
			assert(second->exons[i+1] - second->exons[i] >  0) ;
		}

		fprintf(MY_OUT_FP, "\t0");
		for (uint32_t i = 2; i < second->exons.size(); i+=2)
		{
			fprintf(MY_OUT_FP, ",%d", second->exons[i] - second->exons[0]);
			assert(second->exons[i] - second->exons[0] >  0 ) ;
		}

		fprintf(MY_OUT_FP, "\tqpalmaScore=%1.3f;numMatches=%i;numGaps=%i;minExonLen=%i;maxIntronLen=%i;readOrientation=%c;read=%s",	second->qpalma_score, second->num_matches, second->num_gaps, second->min_exon_len, best->max_intron_len, read_orientation, second->read_anno) ;

		if (_config.POLYTRIM_STRATEGY)
		{
			if (second->polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimStart=%i", second->polytrim_cut_start) ;
			if (second->polytrim_cut_start)
				fprintf(MY_OUT_FP, ";polytrimEnd=%i", second->polytrim_cut_end) ;
		}

    }


	fprintf(MY_OUT_FP, "\n");	
	
	return top_alignments.size() ;
}

int TopAlignments::print_top_alignment_records_shorebed(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP)
{
	int printed = 0 ;
	for (unsigned int i=0; i<top_alignments.size(); i++)
    {
		if (!top_alignments[i]->spliced)
			printed+= print_alignment_shorebed(read, OUT_FP, top_alignments[i], top_alignments[i]->num)  ;
		else
			printed+= print_alignment_shorebed(read, SP_OUT_FP, top_alignments[i], top_alignments[i]->num)  ;
    }
	return printed ;
}

//TODO: merge jörg: check why those shore functions are renamed
int TopAlignments::print_alignment_shorebed(Read const &read, std::ostream *OUT_FP, alignment_t* align, unsigned int num)
{
	HIT* hit=align->hit;
	if (_config.STATISTICS)
		_stats.HITS_MM[hit->mismatches]++;

	int fstart, flen;

	int hitlength = hit->end - hit->start + 1;
	unsigned int readstart;
	if (hit->orientation == '+') {
		readstart = hit->start - hit->readpos + hit->start_offset; // start pos of read in genome	0-initialized
	} else {
		readstart = hit->start - (((int)read.length()) - hit->readpos - hitlength + 2)
			+ hit->start_offset; // 0-initialized
	}

	char FLANK_SEQ[Config::MAX_READ_LENGTH + 200];


	fprintf(OUT_FP, "%s\t%d\t%s\t%s\t%c",
			hit->chromosome->desc(),
			readstart + 1, 	// 1-initialized
			//		ALIGNSEQ,	//Alignment String: ALIGNSEQ is the last recorded alignment string
			align->read_anno,	//Alignment String corresponding to this particular alignment
			read.id(),
			((hit->orientation == '+') ? 'D' : 'P'));

	if (_config.SCORES_OUT)
		fprintf(OUT_FP, "\t%.1f", (double) (hit->gaps * _config.GAP_SCORE
											+ (hit->mismatches - hit->gaps) * _config.MM_SCORE
											- (((int)read.length()) - hit->mismatches) * _config.M_SCORE));
	else
		fprintf(OUT_FP, "\t%d", hit->mismatches);

	fprintf(OUT_FP, "\t%d\t%d\t%d",
			num, 			// Number of hits
			((int)read.length()),	// length of hit on genome
			0);

	if (read.format() == 2)
		fprintf(OUT_FP, "\t%d", read.pe_flag());
	if (strlen(read.quality(0)) != 0)
		fprintf(OUT_FP, "\t%s", read.quality(0));
	if (strlen(read.quality(1)) != 0)
		fprintf(OUT_FP, "\t%s", read.quality(1));
	if (strlen(read.quality(2)) != 0)
		fprintf(OUT_FP, "\t%s", read.quality(2));

	if (_config.FLANKING != 0) {

		fstart = (readstart < _config.FLANKING) ? 0 : readstart - _config.FLANKING;
		flen = (readstart + 
				((int)read.length()) + _config.FLANKING > hit->chromosome->length()) ?
			hit->chromosome->length() - fstart
			: 	readstart - fstart + ((int)read.length())
			+ _config.FLANKING;
		{
			for (size_t i=0; i<(size_t)flen; i++)
				FLANK_SEQ[i]= (*hit->chromosome)[fstart + i];
		}
		FLANK_SEQ[flen] = '\0';
		fprintf(OUT_FP, "\t%d\t%s",
				hit->chromosome->length(),
				FLANK_SEQ);

	} else if (_config.PRINT_SEQ > 0) {
	
		fprintf(OUT_FP, "\t%d", hit->chromosome->length());
	
	}

	fprintf(OUT_FP, "\n");

	return 1;
}

/* In the following (commented) ALIGNSEQ is generated a second time. Since by
   creation of ALIGNSEQ the hit structure is changed (gap encoding), this won't work.
   The function above just uses the already generated ALIGNSEQ variable.
*/
// PERFECT HITS:
/*	if (hit->mismatches == 0) {

	if (hit->orientation == '+')
	strcpy(ALIGNSEQ, _read.data());
	else
	{
	for (size_t i=0; i<_read.length(); i++)
	ALIGNSEQ[i] = (*hit->chromosome)[readstart+i];
	//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + readstart, ((int)_read.lenght()));
	}
	ALIGNSEQ[((int)_read.length())] = '\0';


	// print in file:
	if (_config.OUTPUT_FORMAT == 0) {
	/////// SHORE file ///////

	fprintf(OUT_FP, "%s\t%d\t%s\t%s\t%c", hit->chromosome->desc(),
	readstart + 1, // 1-initialized
	ALIGNSEQ,//Alignment String
	_read.id(), ((hit->orientation == '+') ? 'D' : 'P'));

	if (_config.SCORES_OUT)
	fprintf(OUT_FP, "\t%.1f", (double) (-((int)_read.length()) * _config.M_SCORE));
	else
	fprintf(OUT_FP, "\t%d", ((int)_read.length()) - hit->mismatches);
	// lt - changed this to report num matches instead of num mismatches

	fprintf(OUT_FP, "\t%d\t%d\t%d", num, // Number of hits
	((int)_read.length()), 0);

	if (_read.format() == 2)
	fprintf(OUT_FP, "\t%d", _read.pe_flag());
	if (strlen(_read.quality()[0]) != 0)
	fprintf(OUT_FP, "\t%s", _read.quality()[0]);
	if (strlen(_read.quality()[1]) != 0)
	fprintf(OUT_FP, "\t%s", _read.quality()[1]);
	if (strlen(_read.quality()[2]) != 0)
	fprintf(OUT_FP, "\t%s", _read.quality()[2]);

	if (_config.FLANKING != 0) {
	fstart = (readstart < _config.FLANKING) ? 0 : readstart - _config.FLANKING;
	flen
	= (readstart + ((int)_read.length()) + _config.FLANKING
	> hit->chromosome->length()) ? hit->chromosome->length()
	- fstart
	: readstart - fstart + ((int)_read.length()) + _config.FLANKING;
	{
	for (size_t i=0; i<(size_t)flen; i++)
	FLANK_SEQ[i]=(*hit->chromosome)[fstart + i];
	//strncpy(FLANK_SEQ, CHR_SEQ[hit->chromosome] + fstart, flen);
	}
	FLANK_SEQ[flen] = '\0';
	fprintf(OUT_FP, "\t%d\t%s", hit->chromosome->length(),
	FLANK_SEQ);
	} else if (_config.PRINT_SEQ > 0)
	fprintf(OUT_FP, "\t%d", hit->chromosome->length());
	//if (_config.PRINT_SEQ == 2)
	//	fprintf(OUT_FP, "\t%s", CHR_SEQ[hit->chromosome]);

	fprintf(OUT_FP, "\n");
	} else {
	/////// BED file ///////

	if (_config.SCORES_OUT) {
	fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%.1f\t%c\n",
	hit->chromosome->desc(), readstart, readstart + 1
	+ ((int)_read.length()), _read.id(), (double) (-((int)_read.length())
	* _config.M_SCORE), hit->orientation);
	} else {
	fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%d\t%c\t0\t0\t0,0,0\t0\t0\t0\n",
	hit->chromosome->desc(), readstart, readstart + 1
	+ ((int)_read.length()), _read.id(), ((int)_read.length()) - hit->mismatches,
	hit->orientation);
	// lt  - changed this to report num matches instead of num mismatches
	}

	}

	}
	// HITS WITH MISMATCHES:
	else {

	int count_char = 0;
	char gap_offset = 0;
	char gap_in_read = 0;
	char gap_in_chr = 0;

	// sort mismatches in ascending order according to their abs positions and 'gap before mm'-strategy if equal
	qsort(hit->edit_op, hit->mismatches, sizeof(EDIT_OPS), compare_editops);

	ALIGNSEQ[0] = '\0';

	for (j = 0; j != hit->mismatches; ++j) {

	if (hit->edit_op[j].pos < 0) {
	hit->edit_op[j].pos = -hit->edit_op[j].pos;
	gap_in_chr = 1;
	}

	if (j == 0) 
	{
	for (size_t i=0; (int)i< hit->edit_op[0].pos - 1; i++)
	ALIGNSEQ[i]=(*hit->chromosome)[readstart + i];
	//strncpy(ALIGNSEQ, CHR_SEQ[hit->chromosome] + (readstart),
	//		hit->edit_op[0].pos - 1);
	count_char += hit->edit_op[0].pos - 1;
	} else if (hit->edit_op[j].pos - hit->edit_op[j - 1].pos != 0) 
	{
	for (size_t i=0; (int)i<hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1 + gap_in_read; i++)
	ALIGNSEQ[count_char+i] = (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read) + i];
	//strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome]
	//		+ (readstart + hit->edit_op[j - 1].pos + gap_offset
	//				- gap_in_read), hit->edit_op[j].pos
	//		- hit->edit_op[j - 1].pos - 1 + gap_in_read); // -1???

	count_char += hit->edit_op[j].pos - hit->edit_op[j - 1].pos - 1
	+ gap_in_read;
	} // else: edit_op[j-1] must have been a gap!

	gap_in_read = 0;

	if (hit->edit_op[j].mm) {
	if (hit->orientation == '+')
	sprintf(ALIGNSEQ + count_char, "[%c%c]",
	(*hit->chromosome)[readstart+ hit->edit_op[j].pos - 1 + gap_offset],
	_read.data()[hit->edit_op[j].pos - 1]);
	else
	sprintf(ALIGNSEQ + count_char, "[%c%c]",
	(*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset],
	get_compl_base(_read.data()[((int)_read.length())
	- hit->edit_op[j].pos]));
	} else if (gap_in_chr) {
	if (hit->orientation == '+')
	sprintf(ALIGNSEQ + count_char, "[-%c]",
	_read.data()[hit->edit_op[j].pos - 1]);
	else
	sprintf(ALIGNSEQ + count_char, "[-%c]", get_compl_base(
	_read.data()[((int)_read.length()) - hit->edit_op[j].pos]));

	gap_offset--;
	gap_in_chr = 0;
	} else {
	sprintf(ALIGNSEQ + count_char, "[%c-]",
	(*hit->chromosome)[readstart + hit->edit_op[j].pos - 1 + gap_offset]);

	gap_offset++;
	gap_in_read = 1;
	}

	count_char += 4;

	}

	// from last mismatch to end of read:
	{
	for (size_t i=0; (int)i<((int)_read.length()) - hit->edit_op[j - 1].pos + gap_in_read; i++)
	ALIGNSEQ[count_char+i] = (*hit->chromosome)[(readstart + hit->edit_op[j - 1].pos + gap_offset - gap_in_read) + i];
	//strncpy(ALIGNSEQ + count_char, CHR_SEQ[hit->chromosome] + (readstart
	//	+ hit->edit_op[j - 1].pos + gap_offset - gap_in_read),
	//		  ((int)_read.lenght()) - hit->edit_op[j - 1].pos + gap_in_read);
	}
	count_char += ((int)_read.length()) - hit->edit_op[j - 1].pos + gap_in_read;

	ALIGNSEQ[count_char] = '\0';

	// print in file:
	if (_config.OUTPUT_FORMAT == 0) {

	/////// SHORE file ///////

	fprintf(OUT_FP, "%s\t%d\t%s\t%s\t%c", hit->chromosome->desc(),
	readstart + 1, // 1-initialized
	ALIGNSEQ,//Alignment String
	_read.id(), ((hit->orientation == '+') ? 'D' : 'P'));

	if (_config.SCORES_OUT)
	fprintf(OUT_FP, "\t%.1f", (double) (hit->gaps * _config.GAP_SCORE
	+ (hit->mismatches - hit->gaps) * _config.MM_SCORE
	- (((int)_read.length()) - hit->mismatches) * _config.M_SCORE));
	else
	fprintf(OUT_FP, "\t%d", ((int)_read.length()) - hit->mismatches);
	// lt - changed this to report num matches instead of num mismatches

	fprintf(OUT_FP, "\t%d\t%d\t%d", num, // Number of hits
	((int)_read.length()), // length of hit on genome
	0);

	if (_read.format() == 2)
	fprintf(OUT_FP, "\t%d", _read.pe_flag());
	if (strlen(_read.quality()[0]) != 0)
	fprintf(OUT_FP, "\t%s", _read.quality()[0]);
	if (strlen(_read.quality()[1]) != 0)
	fprintf(OUT_FP, "\t%s", _read.quality()[1]);
	if (strlen(_read.quality()[2]) != 0)
	fprintf(OUT_FP, "\t%s", _read.quality()[2]);

	if (_config.FLANKING != 0) {
	fstart = (readstart < _config.FLANKING) ? 0 : readstart - _config.FLANKING;
	flen
	= (readstart + ((int)_read.length()) + _config.FLANKING
	> hit->chromosome->length()) ? hit->chromosome->length()
	- fstart
	: readstart - fstart + ((int)_read.length()) + _config.FLANKING;
	{
	for (size_t i=0; i<(size_t)flen; i++)
	FLANK_SEQ[i]= (*hit->chromosome)[fstart + i];
	//strncpy(FLANK_SEQ, CHR_SEQ[hit->chromosome] + fstart, flen);
	}
	FLANK_SEQ[flen] = '\0';
	fprintf(OUT_FP, "\t%d\t%s", hit->chromosome->length(),
	FLANK_SEQ);
	} else if (_config.PRINT_SEQ > 0)
	fprintf(OUT_FP, "\t%d", hit->chromosome->length());
	//if (_config.PRINT_SEQ == 2)
	//	fprintf(OUT_FP, "\t%s", CHR_SEQ[hit->chromosome]);

	fprintf(OUT_FP, "\n");
	} else {

	/////// BED file ///////

	if (_config.SCORES_OUT) {
	fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%.1f\t%c\n",
	hit->chromosome->desc(), readstart, readstart + 1
	+ ((int)_read.length()) + gap_offset, _read.id(),
	(double) (hit->gaps * _config.GAP_SCORE + (hit->mismatches
	- hit->gaps) * _config.MM_SCORE - (((int)_read.length())
	- hit->mismatches) * _config.M_SCORE), hit->orientation);
	} else {
	fprintf(OUT_FP, "%s\t%d\t%d\t%s\t%d\t%c\t0\t0\t0,0,0\t0\t0\t0\n",
	hit->chromosome->desc(), readstart, readstart + 1
	+ ((int)_read.length()) + gap_offset, _read.id(),
	((int)_read.length()) - hit->mismatches, hit->orientation);
	// lt - changed this to report num matches instead of num mismatches
	}
	}
	}

	return 1;
	}*/

void TopAlignments::print_bam_header(Genome& genome, FILE*OUT_FP)
{
	for (unsigned int i=0; i<genome.nrChromosomes(); i++)
	{
		fprintf(OUT_FP, "@SQ\tSN:%s\tLN:%i\n", genome.get_desc(i), genome.chromosome(i).length())  ;
		//fprintf(stdout, "@SQ\tSN:%s\tLN:%i\n", genome.get_desc(i), genome.chromosome(i).length())  ;
	}
}

FILE* TopAlignments::open_bam_pipe(std::string & out_fname)
{

	std::string command = std::string("samtools view -Sb /dev/stdin  2> /dev/null > ") + out_fname  + " && echo samtools subprocess terminated successfully";
	if ( _config.OUTPUT_FORMAT_OPTION == OUTPUT_FORMAT_OPTION_SORTPOS )
		command = std::string("samtools view -Sb /dev/stdin 2> /dev/null | samtools sort /dev/stdin ") + out_fname + std::string("&& mv ") + 
			out_fname + ".bam " + out_fname + " && echo samtools subprocess terminated successfully" ;
	if ( _config.OUTPUT_FORMAT_OPTION == OUTPUT_FORMAT_OPTION_SORTNAME )
		command = std::string("samtools view -Sb /dev/stdin  2> /dev/null | samtools sort -n /dev/stdin ") + out_fname + std::string("&& mv ") + 
			out_fname + ".bam " + out_fname + " && echo samtools subprocess terminated successfully" ;
	FILE* OUT_FP=NULL ;
	
	fflush(stdout) ;
	if (false) // does not work
	{
		fprintf(stdout, "Testing samtools pipe: ") ;
		OUT_FP = popen(command.c_str(), "w") ;
		//TopAlignments::print_bam_header(genome, OUT_FP) ;
		int status = pclose(OUT_FP);
		if WIFEXITED(status) 
			{
				int ret=WEXITSTATUS(status);
				if (ret==0)
					fprintf(stdout, "OK\n");
				else
				{
					fprintf(stdout, "exit code %i\n", ret) ;
					fprintf(stderr, "samtools pipe failed (command '%s')\n", command.c_str()) ;
					exit(-1) ;
				}
			}
		fprintf(stderr, "samtools pipe status %i", status) ;
	}
	OUT_FP = popen(command.c_str(), "w") ;
	return OUT_FP ;
}

int TopAlignments::close_bam_pipe(FILE * FP)
{
	int status = pclose(FP);

	if WIFEXITED(status) 
		{
			int ret=WEXITSTATUS(status);
			if (ret!=0)
			{
				return ret ;
			}
			return 0 ;
		}
	return -1 ;
}


// SAM format
int TopAlignments::print_top_alignment_records_sam(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP)
{
	std::ostream* MY_OUT_FP = OUT_FP ;
	if (top_alignments.size()==0){
		
		if (_config.INCLUDE_UNMAPPED_READS_SAM)
		{
			_stats.alignment_num_unmapped+= 1 ;
			_stats.print_alignment_stats() ;

			uint32_t flag=4 ;			
			char unavailablechar='*';
			int  unavailableint=0;
			int mapq=255;
			
			fprintf(MY_OUT_FP, "%s\t%d\t%c\t%i\t%i\t%c\t%c\t%i\t%i", 
					read.id(),
					flag, 
					unavailablechar,
					unavailableint,
					mapq,
					unavailablechar,
					unavailablechar,
					unavailableint,
					unavailableint);
		
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_READ)
				fprintf(MY_OUT_FP, "\t%s", read.data()) ;
			else
				fprintf(MY_OUT_FP, "\t*") ;
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_QUALITY)
				fprintf(MY_OUT_FP, "\t%s", read.quality(0)) ;
			else
				fprintf(MY_OUT_FP, "\t*") ;

			fprintf(MY_OUT_FP, "\n");
			
		}
		return 0 ;
	}
	
	// if (_config.RTRIM_STRATEGY)
	//     return 0;

	// pre compute H0, H1, and H2 tags
	uint32_t H0 = 0 ;
	uint32_t H1 = 0 ;
	uint32_t H2 = 0 ;
	uint32_t min_edit_ops = 1000 ;
	double max_qpalma_score = -1000 ;

	for (unsigned int j=0; j<top_alignments.size(); j++)
    {
		if (top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches == 0)
			H0 += 1 ;
		if (top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches == 1)
			H1 += 1 ;
		if (top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches == 2)
			H2 += 1 ;
		if (top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches < min_edit_ops)
			min_edit_ops = top_alignments[j]->num_gaps + top_alignments[j]->num_mismatches ;
		if (top_alignments[j]->qpalma_score > max_qpalma_score)
			max_qpalma_score = top_alignments[j]->qpalma_score ;
    }

	Read const * curr_read;

	for (unsigned int j=0; j<top_alignments.size(); j++)
    {
		alignment_t * curr_align  = top_alignments[j] ;
		if (curr_align->spliced)
			MY_OUT_FP = SP_OUT_FP ;
		else
			MY_OUT_FP = OUT_FP ;
		
		assert(curr_align->exons.size() >= 2);
		
		int min_exon_len = curr_align->exons[1]-curr_align->exons[0] ;
		int max_intron_len = 0 ;
		int min_intron_len = INT_MAX ;
		for (unsigned int k=2; k<curr_align->exons.size(); k+=2)
		{
			if (min_exon_len>curr_align->exons[k+1]-curr_align->exons[k])
				min_exon_len = curr_align->exons[k+1]-curr_align->exons[k] ;
			if (max_intron_len<curr_align->exons[k]-curr_align->exons[k-1])
				max_intron_len = curr_align->exons[k]-curr_align->exons[k-1] ;
			if (min_intron_len>curr_align->exons[k]-curr_align->exons[k-1])
				min_intron_len = curr_align->exons[k]-curr_align->exons[k-1] ;
		}

		if (j == 0) 
		{
			if (curr_align->rtrim_cut!=0)
				assert(curr_align->polytrim_cut_start==0 && curr_align->polytrim_cut_end==0) ;
			if (curr_align->polytrim_cut_start!=0)
				assert(curr_align->polytrim_cut_end==0) ;
			
			pthread_mutex_lock( &_stats.alignment_num_mutex) ;
			
            if (_config._personality == Palmapper && curr_align->spliced)
            {
                assert(_config.SPLICED_HITS) ;
                _stats.alignment_num_spliced_best++ ;
            } else
                _stats.alignment_num_unspliced_best++ ;
			
			if (_config._personality == Palmapper) {
				for (int i=1; i<(int)top_alignments.size(); i++)
					if (top_alignments[i]->spliced)
						_stats.alignment_num_spliced_suboptimal+=1 ;
					else
						_stats.alignment_num_unspliced_suboptimal+= 1 ;
				_stats.print_alignment_stats() ;
			}
			pthread_mutex_unlock( &_stats.alignment_num_mutex ) ;
		}

		int polytrim_cut_start=curr_align->polytrim_cut_start ;
		int polytrim_cut_end=curr_align->polytrim_cut_end ;
		
		if (_config.RTRIM_STRATEGY && !_config.POLYTRIM_STRATEGY)
		{
			assert(polytrim_cut_end==0) ;
			polytrim_cut_end = curr_align->rtrim_cut ;
		}

		if ((_config.POLYTRIM_STRATEGY || _config.RTRIM_STRATEGY) && (polytrim_cut_start>0 || polytrim_cut_end>0))
			if (read.get_orig())
				curr_read = read.get_orig() ;
			else
				curr_read = &read ;
		else
			curr_read = &read ;

		fprintf(MY_OUT_FP, "%s", curr_align->read_id) ;
		uint32_t flag=0 ;

		// Read orientation output always depends on the left strand -> reverse orientation if found on negative strand
		
		char transcription_direction;
		char read_orientation;
		determine_transcription_direction(curr_align->strand,curr_align->orientation,_config.STRAND,transcription_direction,read_orientation);
		
		//Read reversed and complemented compared to the forward strand
		flag+=((read_orientation=='-')*16) ;

		//For several fragments: First corresponds to left reads and last corresponds to right reads
		flag+=((_config.STRAND==1)*64) ;
		flag+=((_config.STRAND==0)*128) ;

		/* flag+=_config.SEQUENCING_WAS_PAIRED ;
		 * flag+=(curr_read.MAPPED_AS_PAIR*2) ;
		 * flag+=(IS_UNMAPPED*4) ;
		 * flag+=(curr_read.MATE_IS_UNMAPPED*8) ;
		 * flag+=(curr_read.STRAND_OF_MATE*32) ;
		 * flag+=(curr_read.FIRST_IN_PAIR)?64:128 ;
		 */
		flag+=((top_alignments.size()>1)*256) ;

		fprintf(MY_OUT_FP, "\t%d\t%s\t%d", 
				flag, 
				curr_align->chromosome->desc(),
				curr_align->exons[0] + 1) ; 
		if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_MAQQUALITY)
		{
			fprintf(stderr, "MAQ quality not implemented yet\n") ;
			exit(-1) ;
		}
		else
		{
			if (j<=254)
				fprintf(MY_OUT_FP, "\t%i", 254 - j);
			else
				fprintf(MY_OUT_FP, "\t0");
		}
		
		//	double qpalma_score = best->qpalma_score ;
        
		// determine CIGAR
		char __cigar[1000] ; 
		char cigar[1000] ;
		char cig_buf[255] ;

		uint32_t pos_in_cigar = 0;
		uint32_t idx = 0 ;
        
		// handle trimmed start as soft clips
		if ((_config.POLYTRIM_STRATEGY  || _config.RTRIM_STRATEGY) && (polytrim_cut_start>0 && read_orientation=='+'))
		{
			snprintf (cig_buf, (size_t) 255, "%d", polytrim_cut_start) ;
			for (uint32_t ii=0; ii < strlen(cig_buf); ii++)
				cigar[pos_in_cigar + ii] = cig_buf[ii] ;
			pos_in_cigar += strlen(cig_buf) ;
			cigar[pos_in_cigar++] = 'S' ;
		}
		if ((_config.POLYTRIM_STRATEGY  || _config.RTRIM_STRATEGY) && (polytrim_cut_end>0 && read_orientation=='-'))
		{
			snprintf (cig_buf, (size_t) 255, "%d", polytrim_cut_end) ;
			for (uint32_t ii=0; ii < strlen(cig_buf); ii++)
				cigar[pos_in_cigar + ii] = cig_buf[ii] ;
			pos_in_cigar += strlen(cig_buf) ;
			cigar[pos_in_cigar++] = 'S' ;
		}

		//fprintf(stdout,"read anno: %s\n",curr_align->read_anno);
		
		for (uint32_t i = 0; i < strlen(curr_align->read_anno); i++)
		{
			if (curr_align->read_anno[i] != '[')
				__cigar[idx] = 'M' ;	
			else
			{
				if (curr_align->read_anno[i+1] == '-')
					__cigar[idx] = 'I' ;
				else if (curr_align->read_anno[i+2] == '-')
					__cigar[idx] = 'D' ;
				else
					__cigar[idx] = 'M' ;
				i += 3 ;
			}
			idx += 1 ;
		}
		__cigar[idx] = 0 ;


		uint32_t last = __cigar[0] ;
		uint32_t count = 1 ;
		uint32_t ii = 0;
		uint32_t insertions = 0 ;
		uint32_t deletions = 0 ;
		idx = 0 ;
		uint32_t exon_size = (curr_align->exons[idx + 1] - curr_align->exons[idx]) ;


		for (uint32_t i = 1; i < strlen(__cigar); i++)
		{
			// we reached end of current exon, which is not the last exon -> add intron
			if ((i - insertions) == exon_size && idx + 2 < curr_align->exons.size())
			{
				// finish current block
				snprintf (cig_buf, (size_t) 255, "%d", count) ;
				for (ii=0; ii < strlen(cig_buf); ii++)
					cigar[pos_in_cigar + ii] = cig_buf[ii] ;
				pos_in_cigar += strlen(cig_buf) ;
				cigar[pos_in_cigar++] = last ;

				if (last == 'D')
					deletions += count;
				if (last == 'I')
					insertions += count;

				// add intron
				count = 0 ;
				last = ' ' ;
				snprintf (cig_buf, (size_t) 255, "%d", curr_align->exons[idx + 2] - curr_align->exons[idx + 1]) ;
				for (ii=0; ii < strlen(cig_buf); ii++)
					cigar[pos_in_cigar + ii] = cig_buf[ii] ;
				pos_in_cigar += strlen(cig_buf) ;
				cigar[pos_in_cigar++] = 'N' ;

				/// go to next exon
				idx += 2 ;
				exon_size += (curr_align->exons[idx + 1] - curr_align->exons[idx]) ;
				assert(curr_align->exons[idx+1] - curr_align->exons[idx] >  0) ;
			}

			// start new block in cigar
			if (__cigar[i] != (char)last) 
			{
				if (last != ' ')
				{
					snprintf (cig_buf, (size_t) 255, "%d", count) ;
					for(ii=0; ii < strlen(cig_buf); ii++)
						cigar[pos_in_cigar + ii] = cig_buf[ii] ;
					pos_in_cigar += strlen(cig_buf) ;
					cigar[pos_in_cigar++] = last ;
				}
				if (last == 'D')
					deletions += count;
				if (last == 'I')
					insertions += count;
				count = 1 ;
				last = __cigar[i] ;
			}
			else
				count += 1 ;
		}

		// handle last block
		snprintf (cig_buf, (size_t) 255, "%d", count) ;
		for (ii=0; ii < strlen(cig_buf); ii++)
			cigar[pos_in_cigar + ii] = cig_buf[ii] ;
		pos_in_cigar += ii ;
		cigar[pos_in_cigar++] = last ;
		if (last == 'D')
			deletions += count ; 
		if (last == 'I')
			insertions += count ; 


		// handle trimmed reads end
		if ((_config.POLYTRIM_STRATEGY || _config.RTRIM_STRATEGY) && polytrim_cut_end>0 &&  read_orientation=='+')
		{
			snprintf (cig_buf, (size_t) 255, "%d", polytrim_cut_end) ;
			for (ii=0; ii < strlen(cig_buf); ii++)
				cigar[pos_in_cigar + ii] = cig_buf[ii] ;
			pos_in_cigar += strlen(cig_buf) ;
			cigar[pos_in_cigar++] = 'S' ;
		}
		// handle trimmed reads end
		if ((_config.POLYTRIM_STRATEGY || _config.RTRIM_STRATEGY) && polytrim_cut_start>0  && read_orientation=='-')
		{
			snprintf (cig_buf, (size_t) 255, "%d", polytrim_cut_start) ;
			for (ii=0; ii < strlen(cig_buf); ii++)
				cigar[pos_in_cigar + ii] = cig_buf[ii] ;
			pos_in_cigar += strlen(cig_buf) ;
			cigar[pos_in_cigar++] = 'S' ;
		}

		//cigar[pos] = 0 ;
		//if (cum_size + indel_offset + polytrim_cut_start + polytrim_cut_end != curr_read->length()) 
		if (exon_size + insertions - deletions != curr_read->length()-polytrim_cut_start-polytrim_cut_end){
			
			fprintf(stdout, "WARNING - block sum does not match readlength: block_sum=%i, readlength=%i, read=%s, read_id=%s \n", exon_size + insertions - deletions, curr_read->length(), curr_read->data(), curr_align->read_id) ;
			//fprintf(stderr, " cigar: %s\nexon %i, insert %i, deletion %i, trim_start %i, trim_end %i, read_length %i\n read %s\n read anno %s\n\n", cigar, exon_size,insertions,deletions, curr_align->polytrim_cut_start, curr_align->polytrim_cut_end, curr_read->length(), curr_read->data(),curr_align->read_anno) ;
		}
		//assert(cum_size + indel_offset + curr_align->polytrim_cut_start + curr_align->polytrim_cut_end == curr_read->length()) ;

		cigar[pos_in_cigar] = 0 ;


		//		if (curr_align->orientation=='+' || curr_align->exons.size() < 3)
		fprintf(MY_OUT_FP, "\t%s\t*\t0\t0", cigar) ; 
		/*		else
				{
				// reverse order of cigar
				char rcigar[500] ;
				uint32_t marker = 0 ;
				for (int k=1; k<strlen(cigar); k++)
				{
				if (cigar[strlen(cigar)-k-1] <= '9')
				continue ;
				else
				{
				for (uint32_t kk = 0; kk < k - marker; kk++)
				rcigar[marker+kk]=cigar[strlen(cigar)-k+kk] ;
				marker = k;
				}
				}
				for (uint32_t kk = 0; kk < strlen(cigar) - marker; kk++)
				rcigar[marker+kk]=cigar[kk] ;
            
				rcigar[strlen(cigar)]=0 ;
				fprintf(MY_OUT_FP, "\t%s\t*\t0\t0", rcigar) ; 
				}*/

		if (read_orientation=='+')
		{
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_READ)
				fprintf(MY_OUT_FP, "\t%s", curr_read->data()) ;
			else
				fprintf(MY_OUT_FP, "\t*") ;
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_QUALITY)
				fprintf(MY_OUT_FP, "\t%s", curr_read->quality(0)) ;
			else
				fprintf(MY_OUT_FP, "\t*") ;
		}		
		else
		{
			// reverse order of quality
			char qual[500] ;
			for (int k=0; k<((int)curr_read->length()); k++)
				qual[k]=(curr_read->quality(0))[((int)(curr_read->length()))-k-1] ;
			qual[((int)(curr_read->length()))]=0 ;
            
			// complementary reverse read 
			char cr_read[500] ;
			for (int i=0; i<((int)curr_read->length()); i++)
				cr_read[i] = get_compl_base(curr_read->data()[((int)(curr_read->length()))-i-1]) ;
			
			cr_read[((int)(curr_read->length()))]=0 ;
			
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_READ)
				fprintf(MY_OUT_FP, "\t%s", cr_read      );
			else
				fprintf(MY_OUT_FP, "\t*") ;
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_QUALITY)
				fprintf(MY_OUT_FP, "\t%s",qual);
			else
				fprintf(MY_OUT_FP, "\t*") ;
		}

		if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_SAMFLAGS)
			fprintf(MY_OUT_FP, "\tNM:i:%i", curr_align->num_mismatches + curr_align->num_gaps) ;
		if (H0 > 0 && (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_SAMFLAGS))
			fprintf(MY_OUT_FP, "\tH0:i:%i", H0) ;
		if (H1 > 0 && (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_SAMFLAGS))
			fprintf(MY_OUT_FP, "\tH1:i:%i", H1) ;
		if (H2 > 0 && (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_SAMFLAGS))
			fprintf(MY_OUT_FP, "\tH2:i:%i", H2) ;


		if (curr_align->spliced)
		{
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_MORESAMFLAGS)
			{
				// if (!curr_align->non_consensus_alignment)
				// 	fprintf(MY_OUT_FP, "\tXS:A:%c", curr_align->strand) ;
				// else
				// {
				// 	if (_config.STRAND > -1) {
				// 		if ((( curr_align->orientation == '+') && _config.STRAND) || ((curr_align->orientation == '-') && ! _config.STRAND))
				// 			fprintf(MY_OUT_FP, "\tXS:A:+") ;
				// 		else
				// 			fprintf(MY_OUT_FP, "\tXS:A:-") ;
				// 	}
					
				// }
				if(!curr_align->non_consensus_alignment || _config.STRAND > -1) 
					fprintf(MY_OUT_FP, "\tXS:A:%c",transcription_direction) ;

				fprintf(MY_OUT_FP, "\tXe:i:%i", min_exon_len) ;
				fprintf(MY_OUT_FP, "\tXI:i:%i", max_intron_len) ;
				fprintf(MY_OUT_FP, "\tXi:i:%i", min_intron_len) ;
				if (curr_align->non_consensus_alignment) 
				{
					assert(curr_align->intron_consensus.size()>0);
					fprintf(MY_OUT_FP, "\tXC:Z:%s", curr_align->intron_consensus[0].c_str()) ;
					for (unsigned int qq=1; qq<curr_align->intron_consensus.size(); qq++)
						fprintf(MY_OUT_FP, "~~~~%s", curr_align->intron_consensus[qq].c_str()) ;
				}
			}
        }
        else if (_config.STRAND > -1 && _config.PROTOCOL != PROTOCOL_UNSTRANDED) {
			if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_MORESAMFLAGS) {
				fprintf(MY_OUT_FP, "\tXS:A:%c",transcription_direction) ;
				// if ((( curr_align->orientation == '+') && _config.STRAND) || ((curr_align->orientation == '-') && ! _config.STRAND))
				// 	fprintf(MY_OUT_FP, "\tXS:A:+") ;
				// else
				// 	fprintf(MY_OUT_FP, "\tXS:A:-") ;
			}
        }

		if (_config.OUTPUT_FORMAT_FLAGS & OUTPUT_FORMAT_FLAGS_MORESAMFLAGS)
		{
			fprintf(MY_OUT_FP, "\tXQ:i:%i", curr_align->qual_mismatches) ;
			fprintf(MY_OUT_FP, "\tXN:i:%i", (int)curr_align->exons.size()/2) ;
			fprintf(MY_OUT_FP, "\tZS:f:%2.3f", curr_align->qpalma_score) ;
			fprintf(MY_OUT_FP, "\tAS:i:%i", (int)(100*curr_align->qpalma_score)) ;
			fprintf(MY_OUT_FP, "\tHI:i:%i", j) ;
			fprintf(MY_OUT_FP, "\tXD:f:%2.3f", max_qpalma_score-curr_align->qpalma_score) ;
			fprintf(MY_OUT_FP, "\tXd:i:%i", curr_align->num_mismatches + curr_align->num_gaps - min_edit_ops) ;
			//fprintf(MY_OUT_FP, "\tXr:i:%i", curr_align->remapped + (curr_align->exons.size()>2)) ;
		}
		fprintf(MY_OUT_FP, "\n") ;

    }

	return top_alignments.size() ;
}

