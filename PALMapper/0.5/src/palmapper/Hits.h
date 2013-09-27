#pragma once

#include <assert.h>
#include <vector>

#include <palmapper/Config.h>
#include <palmapper/Chromosome.h>
#include <palmapper/QPalma.h>
#include <palmapper/QueryFile.h>
#include <palmapper/Read.h>
#include <palmapper/TopAlignments.h>

#define CONTAINER_SIZE 100000
#define SCORE_INTERVAL 1

extern Config _config;

extern char *get_seq(Read const &read, unsigned int n);

typedef struct edit_op_structure {
	signed int pos;
	int mm: 1;
} EDIT_OPS;

struct HIT {
	unsigned short int readpos;
	unsigned int start;
	unsigned int end;
	Chromosome const *chromosome;
	char orientation;
	unsigned char mismatches;
	unsigned char gaps;
	signed char start_offset;
	signed char end_offset;
	EDIT_OPS edit_op[Config::MAX_EDIT_OPS];
	char conversion;	// only for BSSEQ mode
	char aligned;
	HIT *same_eo_succ;	// the list of HITS_BY_SCORE - only forward pointer for now
	HIT *next;
	HIT *last;
};

extern void printhit(Read const &read, HIT* hit);

template <bool clearNew, class Entry> class Container {
public:
	Container() : _entries(30) {
		_containerSize = 64;
		_current = _end = NULL;
		_used = 0;
	}
	~Container() {
		clear();
	}
	Entry *newEntry() {
		if (_current >= _end) {
			_containerSize *= 2;
			_entries.push_back(_current = new Entry[_containerSize]);
			_end = _current + _containerSize;
		}
		if (clearNew)
			::memset(_current, 0, sizeof(Entry));
		++_used;
		return _current++;
	}
	void clear() {
		while (_entries.size() > 0) {
			delete[] _entries.back();
			_entries.pop_back();
		}
		_current = _end = NULL;
		_used = 0;
		_containerSize = 64;
	}

	unsigned int used() {
		return _used;
	}

private:
	std::vector<Entry*> _entries;
	Entry *_current;
	Entry *_end;
	unsigned int _used;
	unsigned int _containerSize;
};

struct MAPPING_ENTRY {
	unsigned int readpos;
	HIT *hit;
	struct MAPPING_ENTRY *pred;
	struct MAPPING_ENTRY *succ;
};

struct CHROMOSOME_ENTRY {
	int chromosome; // @TODO? Minus values indicate reverse hits (would save strand var = 1 byte, ~15MB)
	unsigned int genome_pos;
	char strand;

	CHROMOSOME_ENTRY *next;
	MAPPING_ENTRY *mapping_entries;

	// It seems to be cheaper to store the back-pointer information (pos)
	// in each of these entries rather than having a superior structure.
};


template <class T> class  MyArr {
public:
	MyArr(int size) {
		_accessCount = _accessCountConst = 0;
		_array = new T[size];
		memset(_array, 0, sizeof(T) * size);
	}

	~MyArr() {
		delete[] _array;
	}

	T const &operator[] (size_t index) const {
		return _array[index];
	}

	T &operator[] (size_t index) {
		_used.push_back(index);
		return _array[index];
	}

	T const *operator + (int size) const {
		return &(*this)[size];
	}

	int accesssCountConst() const {
		return _accessCountConst;
	}

	int accesssCount() const {
		return _accessCount;
	}

	void clear() {
		for (size_t i = 0; i < _used.size(); ++i)
			_array[_used[i]] = NULL;
		_used.clear();
	}

private:
	std::vector<size_t> _used;
	mutable int _accessCountConst;
	mutable int _accessCount;
	T *_array;
};

typedef MyArr<CHROMOSOME_ENTRY*> GenomeArr;
//typedef GenomeMap<CHROMOSOME_ENTRY*> GenomeArr;

typedef Container<false, CHROMOSOME_ENTRY> CHROMOSOME_ENTRY_CONTAINER;

typedef struct hits_by_score_structure {
	HIT *hitpointer;
	int num;
} HITS_BY_SCORE_STRUCT;


class Genome ;
class TopAlignments ;
class QPalma ;
class GenomeMaps ;

enum index_type_t { 
	array=1, bwt=2, debug=3
}  ;

	
class Hits {
	friend class Mapper;
public:

	Hits(Genome const &genome, GenomeMaps &genomeMaps, Mapper &hits, Read const &read);

	~Hits();
	Read const &getRead() const {
		return _read;
	}
	int dealloc_hit_lists_operator() ;
	int dealloc_hits_by_score() ;
	int align_hit_simple(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, int orientation, unsigned char mismatches) ;
	int prepare_kbound_alignment(HIT* hit, int start, int end, int readpos, Chromosome const &chromosome, char orientation, char mismatches) ;
	int size_hit(HIT *hit, unsigned int oldlength);
	template<enum index_type_t> int seed2genome(unsigned int num, unsigned int readpos, char conversion);
	void printgenome() ;
//	int alloc_genome_memory() ;
	int duplicate(HIT* hit) ;
	int insert_into_scorelist(HIT* hit, char d) ;
	int browse_hits() ;
	template<enum index_type_t index_type> int map_fast(Read & read);
	template<enum index_type_t index_type> int map_short_read(Read & read, unsigned int num);
	template<enum index_type_t> int get_slots(Read & read, int pos);

	void printhits() ;
	int init_from_meta_index() ;
	int init_hit_lists()  ;

	int report_read_alignment(HIT* hit, int nbest)  ;

	HIT **HIT_LISTS_OPERATOR;
	int analyze_hits(QPalma const * qpalma);
	unsigned int HITS_IN_SCORE_LIST;
	void dealloc_mapping_entries() {_mappings.clear();}
	void dealloc_hits() {_hits.clear();}

	TopAlignments &topAlignments() {
		return _topAlignments;
	}
	int map_reads(Genome &genome, TopAlignments* topalignments) ;

	int get_num_edit_ops() const {
		return _numEditOps;
	}

	inline void reset_num_edit_ops()
	{
		SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.resize(0) ;
	}

	inline void update_num_edit_ops(int num_edit_ops, char & all_hit_strategy, int & NUM_EDIT_OPS_)
	{
		assert(num_edit_ops<=Config::MAX_EDIT_OPS) ;

		//if (!all_hit_strategy && num_edit_ops < NUM_EDIT_OPS_)
		//	NUM_EDIT_OPS_ = num_edit_ops ;

		if (_config.OUTPUT_FILTER==OUTPUT_FILTER_TOP)
		{
			// keep a list of minimal edit operation (assuming that each hit is only reported once)
			// the list is 10 times as long as NUM_TOP_ALIGNMENTS, as will be filtered later according to
			// the qpalma-alignment score (this is a heuristic, which may work well enough; alternatively, one
			// would need to compute the qpalma score for each hit, which will be too expensive

			bool inserted = false ;

			std::vector<int>::iterator it = SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.begin();

			for (uint8_t i = 0; i < SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size(); i++, it++)
			{
				if ( num_edit_ops < SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[i] )
				{
					SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.insert(it, num_edit_ops);
					inserted = true;

					if (SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size() > _config.OUTPUT_FILTER_NUM_TOP*10)
						SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.pop_back();

					break;
				}
			}
			if (!inserted && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size() < _config.OUTPUT_FILTER_NUM_TOP*10)
			{
				SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.push_back(num_edit_ops) ;
				inserted = true ;
			}
		}

		if (!all_hit_strategy && SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()> _config.OUTPUT_FILTER_NUM_TOP*8)
			NUM_EDIT_OPS_ = SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS.size()-1] ;

		//fprintf(stdout, "update_num_edit_ops(num_edit_ops=%i, all_hit_strategy=%i, NUM_EDIT_OPS_=%i)\n", num_edit_ops, (int)all_hit_strategy, NUM_EDIT_OPS_) ;

	//void make_slots(Read read, std::string slotseq, int conversion);
		
	}

	void clear() {
		dealloc_mapping_entries();
		dealloc_hits();
		dealloc_hits_by_score();
		CHROMOSOME_ENTRY_OPERATOR.clear();
		GENOME.clear();
		dealloc_hit_lists_operator();
	}

private:
	int alloc_hits_by_score() ;
	int alloc_hit_lists_operator() ;
	char* get_seq(unsigned int n);
	CHROMOSOME_ENTRY* alloc_chromosome_entry(Read const &read, unsigned int pos, Chromosome const &chr, char strand);
	unsigned int extend_seed(int direction, unsigned int seed_depth_extra, Chromosome const &chr, int genome_pos, int readpos, char conversion);
	int features_conversion(HIT* hit) ;
	int map_fast_bsseq(Read& read, int run, int nr_runs, char conversion);
	void generate_all_possible_seeds(Read & read, int num, int seedpos, unsigned int iter, int fwd_slot, int rev_slot, char conversion);

	Genome const &_genome;
	GenomeMaps &_genomeMaps;
	HITS_BY_SCORE_STRUCT *HITS_BY_SCORE;
	unsigned int NUM_SCORE_INTERVALS;
	Container<true, MAPPING_ENTRY> _mappings;
	Container<true, HIT> _hits;
	Mapper &_mapper;
	Read const &_read;
	GenomeArr &GENOME;
	CHROMOSOME_ENTRY_CONTAINER &CHROMOSOME_ENTRY_OPERATOR;
	std::vector<int> SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS ;
	TopAlignments _topAlignments;
	int _numEditOps;
	char ALL_HIT_STRATEGY;
	unsigned int LONGEST_HIT;
	clock_t seed_covered_reporting_time;

	char HAS_SLOT;
	unsigned int SLOTS[2];
	char *SLOT_STR[2] ;
	int SLOT_STR_POS[2] ;
	std::vector<unsigned int> SLOTS_CV1_FWD;
	std::vector<unsigned int> SLOTS_CV1_REV;
	std::vector<unsigned int> SLOTS_CV2_FWD;
	std::vector<unsigned int> SLOTS_CV2_REV;
};
