#pragma once

#include <palmapper/Config.h>
#include <palmapper/Hits.h>
#include <palmapper/Read.h>
#include <palmapper/dyn_prog/qpalma_dp.h>
#include <palmapper/TopAlignments.h>
#include <palmapper/JunctionMap.h>

struct region_t {
	int32_t start;
	int32_t end;
	bool erased ;
	bool from_map ;
	char orientation ;
	bool* read_map ;
	//int32_t chromosome ;
	//char strand ;
};

class GenomeMaps;
struct HIT;
class QPalma ;
class Hits;
class TopAlignments;

struct alignment_parameter_struct {
	struct penalty_struct h, a, d, *qualityPlifs;
	int num_qualityPlifs;
	double *matchmatrix;
	int matchmatrix_dim[2];
	int quality_offset ;
} ;


class QPalma
{
public:
	class Result {
	public:
		Result(Read const &read, QPalma const &qpalma);
		void add_buffer_to_region(int ori, Chromosome const &chrN, int32_t nregion) ;
		void delete_regions() ;
		void cleanup() {
		  for (int i = 0; i < 2; i++)
			for (int32_t chrN = 0; chrN < (int)regions[i].size(); chrN++) {
			  regions[i][chrN].clear();
			}
		}

	public:
		std::vector<std::vector<region_t *> > regions[2];
		std::vector<std::vector<region_t *> > long_regions[2];
		int qpalma_filter_reason;
		QPalma const &_qpalma;
		Read const &_read;
	};

	struct perform_alignment_t
	{
		QPalma::Result *result;
		Hits *readMappings;
		std::string read_string ;
		std::string read_quality ;
		std::string dna ;
		std::vector<region_t *> current_regions ;
		std::vector<int> positions ;
		Chromosome const *contig_idx ;
		char strand ;
		int ori ;
		int ret ;
		int hit_read;
		int hit_dna;
		int hit_length;
		QPalma const * qpalma ;
		bool joined ;
		bool non_consensus_search ;
		bool remapping;
		JunctionMap *annotatedjunctions;
		ALIGNMENT* aln;
	} ;

    // initialization
	////////////////////

	QPalma(Genome* genome_, GenomeMaps* genomemaps_, int verbosity_=2) ;
	~QPalma() ;

	int map_splice_sites(std::string file_template, char type, float &splice_site_threshold, bool estimate_thresh, bool do_report) ;
protected:
	int init_spliced_align(const char *fname, struct penalty_struct &h,
						   struct penalty_struct &a, struct penalty_struct &d,
						   struct penalty_struct *&qualityPlifs, int &num_qualityPlifs,
						   double*&matchmatrix, int dims[2], int &quality_offset) ;
	int check_splice_files(std::string file_template) ;
	void skip_comment_lines(FILE* fd) ;
	int read_matrix(FILE* fd, double *& matrix, char*& name, int dims[2]) ;
	int read_plif(FILE *fd, struct penalty_struct &plif) ;
	int clean_alignment_parameters() ;
	int init_alignment_parameters(std::string qpalma_file) ;
	int compare_double(const void *a, const void *b) ;
	

    // qpalma filtering
	////////////////////
	
public:
	int qpalma_filter(Result &result, struct alignment_t *ali, int num_N,JunctionMap &annotatedjunctions) const;
	//void qpalma_filter_stat(Result &result, bool spliced) const;
	//void qpalma_filter_stat_report() const;

	int get_qpalma_quality_offset() 
	{
		if (!alignment_parameters)
			return 33 ;
		return alignment_parameters->quality_offset ;
	}
	
	
protected:
	int get_num_splicesites(std::string file_template, const char* type, Chromosome const &chr, char strand, int start, int end, float thresh) const;
	

 public:
	
    // qpalma alignment
	////////////////////

public:
	int capture_hits(Hits &hits, Result &result, bool const non_consensus_search, JunctionMap &annotatedjunctions) const;
	int capture_hits_2(Hits &hits, Result &result, bool const non_consensus_search, JunctionMap &annotatedjunctions) const;
	int junctions_remapping(Hits &hits, Result &result, JunctionMap &junctionmap, int nb_spliced_alignments, JunctionMap &annotatedjunctions) const;
	int perform_alignment(Result &result, Hits &readMappings, std::string &read_string, std::string &read_quality, std::string &dna, std::vector<region_t *> &regions, std::vector<int> &positions,
						  Chromosome const &contig_id, char strand, int ori, int hit_read, int hit_dna, int hit_length, bool non_consensus_search, ALIGNMENT *& aln, bool remapping, JunctionMap &annotatedjunctions) const;
	double score_unspliced(Read const &read, const char * read_anno, const char strand, const char ori) const;
	//void capture_hits_timing(int read_count=-1, float this_read=-1.0) const;
	
protected:
	
	int get_splicesite_positions(std::string file_template, const char *type, Chromosome const &chr, char strand, int start, int end, float thresh, bool store_pos,
								 std::vector<int> &positions) const;
	
	
	int get_string_from_region(Chromosome const &chrN, region_t *region, std::string &str) const;
	void qsort(region_t** output, int size) const;
	void recover_long_regions(Read const &read, std::vector<region_t*> &long_regions_output, std::vector<region_t*> long_regions, std::vector<region_t*> current_regions) const;
	int convert_dna_position(int real_position, size_t* cum_length, const std::vector<region_t *> &current_regions) const;
	int get_first_read_map(Read const &read, bool* read_map) const;
	void print_hit(HIT *hit) ;
	void print_region(region_t *region, const char * bla)  ;
	void print_map(Read const &read, bool* read_map, const char *name) ;

	int perform_alignment_starter(Result &result, Hits &readMappings, std::string read_string, std::string read_quality, std::string dna, std::vector<region_t *> current_regions, std::vector<int> positions, Chromosome const &contig_idx, char strand, int ori, int hit_read_position, int hit_dna_position, int hit_length, bool non_consensus_search, int &num_alignments_reported, bool remapping, JunctionMap &annotatedjunctions) const;

	void delete_long_regions(std::vector<std::vector<region_t *> > *long_regions) const;
	int get_transcription_direction(int side,int orientation) const;
	
//	int rescue_alignment(Read const &read, std::string & read_anno, int ori, int &num_A, int &num_T, int &num) ;


	// inline helpers
	////////////////////

	static inline int ori_map(char c)
	{
		if (c == '+')
			return 0;
		if (c == '-')
			return 1;
		fprintf(stdout, "ori: %c\n", c);
		
		assert(0);
	};

	static inline std::string reverse(std::string str)
	{
		for (int i = 0; i < (int)str.length() / 2; i++) 
		{
			char c = str[i];
			str[i] = str[str.length() - i - 1];
			str[str.length() - i - 1] = c;
		}
		
		return str;
	}
	
	static inline std::vector<int> reverse(std::vector<int> vec)
	{
		for (int i = 0; i < (int)vec.size() / 2; i++) 
		{
			int c = vec[i];
			vec[i] = vec[vec.size() - i - 1];
			vec[vec.size() - i - 1] = c;
		}
	
		return vec;
	}
	
	static inline void reverse(double *vec, int len)
	{
		for (int i = 0; i < (int)len / 2; i++) {
			double c = vec[i];
			vec[i] = vec[len - i - 1];
			vec[len - i - 1] = c;
		}
	}
	
	static inline char complement(char c) 
	{
	  switch (c) {
	  case 'a':
	    return 't';
	  case 'c':
	    return 'g';
	  case 'g':
	    return 'c';
	  case 't':
	    return 'a';
	  case 'A':
	    return 'T';
	  case 'C':
	    return 'G';
	  case 'G':
	    return 'C';
	  case 'T':
	    return 'A';
	  case '[':
	    return ']';
	  case ']':
	    return '[';
	  case '-':
	    return '-';
	  default:
	    if (c >= 'a' && c <= 'z')
	      return 'n';
	    else if (c >= 'A' && c <= 'Z')
	      return 'N';
	    else
	      assert(0);
	  }
	}

	static inline std::string complement(std::string str) {
		for (int i = 0; i < (int)str.length(); i++) {
			char c = str[i];
			switch (c) {
			case 'a':
				str[i] = 't';
				break;
			case 'c':
				str[i] = 'g';
				break;
			case 'g':
				str[i] = 'c';
				break;
			case 't':
				str[i] = 'a';
				break;
			case 'A':
				str[i] = 'T';
				break;
			case 'C':
				str[i] = 'G';
				break;
			case 'G':
				str[i] = 'C';
				break;
			case 'T':
				str[i] = 'A';
				break;
			case '[':
				str[i] = ']';
				break;
			case ']':
				str[i] = '[';
				break;
			case '-':
				str[i] = '-';
				break;
			default:
				if (c >= 'a' && c <= 'z')
					str[i] = 'n';
				else if (c >= 'A' && c <= 'Z')
					str[i] = 'N';
				else
					assert(0);
			}
		}
		
		return str;
	}

protected:
	
	struct alignment_parameter_struct *alignment_parameters;

	const int verbosity ;
	const int MIN_NUM_MATCHES ;
	
	Genome * genome ;
	GenomeMaps* genomemaps ;
} ;
