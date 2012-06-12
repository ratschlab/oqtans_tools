#pragma once

// ##############################################################
// ####### GenomeMaps ###########################################
// ##############################################################

//#define CHR_MAP_DNAARRAY
//#define CHR_MAP_DNAARRAY_CLASS CDNAArray4
// define this if CDNAArray4 is chosen
//#define CHR_MAP_DNAARRAY_2BIT

#include <palmapper/Chromosome.h>

#ifndef CHR_MAP_DNAARRAY
#define MASK_MAPPED_READ_BEST      1
#define MASK_MAPPED_READ           2
#define MASK_SPLICED_READ_BEST      4
#define MASK_SPLICED_READ           8
#else
#ifdef CHR_MAP_DNAARRAY_2BIT
#define MASK_MAPPED_READ_BEST      1
#define MASK_MAPPED_READ           4
#define MASK_SPLICED_READ_BEST      2
#define MASK_SPLICED_READ_BEST_ORIG      4
#define MASK_SPLICED_READ           8
#else
#define MASK_MAPPED_READ_BEST      1
#define MASK_MAPPED_READ           2
#define MASK_SPLICED_READ_BEST      4
#define MASK_SPLICED_READ           8
#endif
#endif
#define MASK_MAPPED_REGION         16
#define MASK_REPETITIVE_SEED       32
#define MASK_REPETITIVE_SEED_MANY1  64
#define MASK_REPETITIVE_SEED_MANY2  128


class GenomeMaps
{
public:
	GenomeMaps(Genome const &genome) ;
	~GenomeMaps() ;

	inline unsigned char CHR_MAP(Chromosome const &chr, size_t index)
	{
#ifdef CHR_MAP_DNAARRAY
		return CHR_MAP_a[chr.nr()]->get_elem(index) ;
#else
		assert(CHR_MAP_c!=NULL) ;
		return CHR_MAP_c[chr.nr()][index] ;
#endif
	}
	

	inline unsigned char CHR_MAP_cov(Chromosome const &chr, size_t index)
	{
#ifdef CHR_MAP_DNAARRAY
		assert(false) ;
		// not implemented
		return CHR_MAP_a[chr.nr()]->get_elem(index) ;
#else
		assert(CHR_MAP_i!=NULL) ;
		return CHR_MAP_i[chr.nr()][index] ;
#endif
	}
	

	inline void CHR_MAP_set(Chromosome const &chr, size_t index, unsigned char c)
	{
#ifdef CHR_MAP_DNAARRAY
#ifdef CHR_MAP_DNAARRAY_2BIT
		assert(c<4) ;
#else // CHR_MAP_DNAARRAY_2BIT
		assert(c<16) ;
#endif // CHR_MAP_DNAARRAY_2BIT
		CHR_MAP_a[chr.nr()]->set_elem(index, c) ;
#else // CHR_MAP_DNAARRAY
		CHR_MAP_c[chr.nr()][index]=c ;
		//CHR_MAP_a[chr]->set_elem(index, c) ;
#endif // CHR_MAP_DNAARRAY
	}

	inline void CHR_MAP_set_cov(Chromosome const &chr, size_t index, unsigned int c)
	{
		CHR_MAP_i[chr.nr()][index]=c ;
	}
	

	int init_reporting() ;
	int report_repetitive_seed(Chromosome const &chr, int chr_start, int count) ;
	int report_mapped_region(Chromosome const &chr, int chr_start, int chr_end, int num_matches)  ;
	int report_mapped_read(Chromosome const &chr, int start, int end, int num_matches, int nbest_hit) ;
	int report_spliced_read(Chromosome const &chr, std::vector<int> & exons, int num_matches, int nbest_hit) ;
	int do_reporting(int force=0) ;
	int read_reporting() ;
	int write_reporting() ;
	int write_cov_reporting() ;
	int clean_reporting() ;
	int init_with_gff(std::string &gff_fname) ;
	
protected:
	unsigned char **CHR_MAP_c ;
	unsigned int **CHR_MAP_i ;
#ifdef CHR_MAP_DNAARRAY
	std::vector<CHR_MAP_DNAARRAY_CLASS*> CHR_MAP_a  ;
#endif
	void to_dnaarray(int chr=-1) ;
	void from_dnaarray(int chr = -1) ;

	// stats
	int reported_repetitive_seeds  ;
	int reported_mapped_regions  ;
	int reported_mapped_reads  ;
	int reported_spliced_reads  ;
	
	int covered_mapped_read_positions  ;
	int covered_mapped_read_positions_best  ;
	int covered_spliced_read_positions  ;
	int covered_spliced_read_positions_best  ;
	int covered_repetitive_seed_positions  ;
	int covered_repetitive_seed_positions_many1  ;
	int covered_repetitive_seed_positions_many2  ;
	int covered_mapped_region_positions  ;

	static clock_t last_report ;

	Genome const *genome;
	
public:
	//TODO: dd check for multithreading
	int REPORT_REPETITIVE_SEED_DEPTH_EXTRA ;

} ;

