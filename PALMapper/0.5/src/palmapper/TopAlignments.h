#pragma once

#include <pthread.h>
#include <sstream>
#include <palmapper/Config.h>
#include <palmapper/GenomeMaps.h>
#include <palmapper/JunctionMap.h>

class Mapper ;
class GenomeMaps ;
class QPalma ;
struct HIT;

typedef struct alignment_t {
	double qpalma_score;
	double sort_key;
	
	uint32_t num_matches;
	uint32_t num_mismatches;
	uint32_t qual_mismatches;
	uint32_t num_gaps;
	char read_anno[4*Config::MAX_READ_LENGTH] ;
	std::vector<int> exons;
	Chromosome const *chromosome;
	char orientation;
	char strand ;
	char read_id[Config::MAX_READ_ID_LENGTH];
	int min_exon_len ;
	int max_intron_len ;
	bool spliced ;
	//bool rescued ;
	//int rescue_start, rescue_end ;
	int rtrim_cut ;
	int polytrim_cut_start ;
	int polytrim_cut_end ;
	char from_gm ;
	bool passed_filters;
	bool remapped ;
	
	
	HIT* hit ;
	int num ;
	std::vector<std::string> intron_consensus ;
	std::vector<bool> non_consensus_intron ;
	bool non_consensus_alignment ;
	
} ALIGNMENT;

class TopAlignments {

public:
    TopAlignments(GenomeMaps* genomemaps_) ;
	~TopAlignments()
	{
		free(ALIGNSEQ);
		clean_top_alignment_record() ;
	}

	u_int8_t report_unspliced_hit(Read const &read, HIT *hit, int num, QPalma const * qpalma)  ;
	int construct_aligned_string(Read const &read, HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *qual_mismatches_p, int *num_matches_p);
	alignment_t *gen_alignment_from_hit(Read const &read, HIT *best_hit, QPalma const *qpalma) ;
	int construct_aligned_string(HIT *hit, int *num_gaps_p, int *num_mismatches_p, int *num_matches_p);
	

	void init_top_alignment_indice()  ;
	void update_top_alignment_indice()  ;
	alignment_t * add_alignment_record(alignment_t *alignment, int num_alignments);
	

	void clean_top_alignment_record()  ;
	void start_top_alignment_record()  ;
	void check_alignment(struct alignment_t * alignment) ;
	void end_top_alignment_record(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP, int rtrim_cut, int polytrim_start, int polytrim_end, JunctionMap &junctionmap) ;

	int print_top_alignment_records(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP) ;
	int print_top_alignment_records_bedx(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP);
	int print_top_alignment_records_shorebed(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP) ;
	int print_alignment_shorebed(Read const &read, std::ostream *OUT_FP, alignment_t* align, unsigned int num)  ;
	int print_top_alignment_records_sam(Read const &read, std::ostream *OUT_FP, std::ostream *SP_OUT_FP) ;

	static void print_bam_header(Genome& genome, FILE *OUT_FP);
	static FILE* open_bam_pipe(std::string & out_fname) ;
	static int close_bam_pipe(FILE * FP) ;

	size_t size()
	{
		return top_alignments.size() ;
	}
	alignment_t * get_alignment(int idx)
	{
		return top_alignments[idx] ;
	}
	

protected:

	int32_t compare_score(alignment_t *a1, alignment_t *a2) ;
	bool alignment_is_equal(alignment_t *a1, alignment_t *a2) ;
	int alignment_is_opposite(alignment_t *a1, alignment_t *a2) ;

	int spliced_is_overlapping(alignment_t *a1, alignment_t *a2) ;
	int non_consensus_overlaps_consensus(alignment_t *a1, alignment_t *a2) ;
	void sort_top_alignment_list()  ;  
	void qsort_top_alignments(alignment_t** output, int size);
	
	void determine_transcription_direction(char strand,char orientation, int side, char &transcription, char &read_forward);
	
	std::vector<alignment_t *> top_alignments;
	int num_spliced_alignments;
	int num_unspliced_alignments;
	const int verbosity ;
	char *ALIGNSEQ;

	int num_filtered;
	int current_ind;
	int temp_ind;

	//pthread_mutex_t top_mutex;

	GenomeMaps* genomemaps ;
} ;
