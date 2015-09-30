#pragma once

#include <string.h>

#include <string>

#include <palmapper/Config.h>
#include <palmapper/QueryFile.h>
#include <palmapper/Statistics.h>
#include <palmapper/Util.h>

class Read {
	friend class QueryFile;
public:
	Read(QueryFile &queryFile);
	Read(Read const &src, unsigned cutStart = 0, unsigned cutEnd = 0);
	~Read();

	unsigned int length() const {
		return READ_LENGTH;
	}

	char *data() {
		return READ;
	}

	char const *data() const {
		return READ;
	}

	char const *id() const {
		return READ_ID;
	}

	char const *quality(int index) const {
		return READ_QUALITY[index];
	}

	int pe_flag() const {
		return READ_PE_FLAG;
	}

	int pe_type() const {
		if (READ_PE_FLAG < 3) return READ_PE_FLAG;
		else if ((READ_PE_FLAG - 3) % 6 < 3) return 1;
		else return 2;
	}

	char format() const {
		return READ_FORMAT;
	}

	void printOn(std::ostream &out) const;

	void cutOffLast() {
		--READ_LENGTH;
		READ[READ_LENGTH] = '\0';
		READ_QUALITY[0][READ_LENGTH] = '\0';
	}

	void cutOffFirst() {
		for (unsigned int i=0; i<READ_LENGTH; i++)
		{
			READ[i] = READ[i+1];
			READ_QUALITY[0][i] = READ_QUALITY[0][i+1];
		}
		--READ_LENGTH;
	}

	void trim(unsigned cutStart, unsigned cutEnd) {
		copyFrom(*this, cutStart, cutEnd);
	}

//	void trim_read_start(Read * read, int trim_start)
//	{
//		for (unsigned int i=trim_start; i<read->length(); i++)
//		{
//			READ[i-trim_start] = read->READ[i];
//			READ_QUALITY[0][i-trim_start] = read->READ_QUALITY[0][i];
//		}
//		READ_LENGTH = read->length()-trim_start ;
//
//		READ[READ_LENGTH] = '\0';
//		READ_QUALITY[0][READ_LENGTH] = '\0';
//	}
//
//	void trim_read_end(Read * read, int trim_end)
//	{
//		for (unsigned int i=0; i<read->length()-trim_end; i++)
//		{
//			READ[i] = read->READ[i];
//			READ_QUALITY[0][i] = read->READ_QUALITY[0][i];
//		}
//		READ_LENGTH = read->length()-trim_end ;
//
//		READ[READ_LENGTH] = '\0';
//		READ_QUALITY[0][READ_LENGTH] = '\0';
//	}

	bool is_match(char a, char b)
		{
			switch(a)
			{
			case 'a':
			case 'A':
				return b=='A' || b=='a' || b=='N' || b=='n' ;
			case 'c':
			case 'C':
				return b=='C' || b=='c' || b=='N' || b=='n' ;
			case 'g':
			case 'G':
				return b=='G' || b=='g' || b=='N' || b=='n' ;
			case 't':
			case 'T':
				return b=='T' || b=='t' || b=='N' || b=='n' ;
			default:
				return toupper(a)==toupper(b) || toupper(a)=='N' || toupper(b)=='N' ;
			}
		}

	int count_matches_fwd(int read_pos, const char* adapter, const int adapter_len, int & len)
		{
			int matches= 0 ;
			for (int i=0; i<read_pos && adapter_len-i-1>=0; i++, len++)
				matches += is_match(READ[read_pos-i-1], adapter[adapter_len-i-1]) ;
			return matches ;
		}

	int count_matches_rev(int read_pos, const char* adapter, const int adapter_len, int & len)
		{
			int matches= 0 ;
			for (int i=0; read_pos+i<(int)READ_LENGTH && i<adapter_len; i++, len++)
				matches += is_match(READ[read_pos+i], adapter[i]) ;
			return matches ;
		}
		
	void find_adapter(unsigned int &adapter_length_start, unsigned int &adapter_length_end, float frac=0.7)
	{
		const char *ad_fwd="ACACTCTTTCCCTACACGACGCTCTTCCGATCT" ;
		const int ad_len_fwd = strlen(ad_fwd) ;
		const char *ad_rev="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG" ;
		const int ad_len_rev = strlen(ad_rev) ;
		
		{
			int best_m=0 ;
			int best_i=0 ;
			int best_len = 0 ;

			for (unsigned int i=1; i<READ_LENGTH/2; i++)
			{
				int len=0 ;
				int m=count_matches_fwd(i, ad_fwd, ad_len_fwd, len) ;
				if (m>best_m)
				{
					best_m=m ;
					best_i=i ;
					best_len=len ;
				}
			}
			
			if (best_len>=6 && (float)best_m/(float)best_len >= frac)
			{
				adapter_length_start = best_i ;
				//fprintf(stdout, "found adapter fwd: read_id=%s\nread=%s\nbest_m=%i\nbest_i=%i\nbest_len=%i\n", READ_ID, READ, best_m, best_i, best_len) ;
				//exit(-1) ;
			}
		}
		
		{
			int best_m=0 ;
			int best_i=0 ;
			int best_len=0 ;
			for (unsigned int i=READ_LENGTH/2; i<READ_LENGTH; i++)
			{
				int len=0 ;
				int m=count_matches_rev(i, ad_rev, ad_len_rev, len) ;
				if (m>best_m)
				{
					best_m=m ;
					best_i=i ;
					best_len=len ;
				}
			}
			
			if (best_len >=6 && (float)best_m/(float)best_len >= frac)
			{
				adapter_length_end = READ_LENGTH-best_i ;
				//fprintf(stdout, "found adapter rev: read_id=%s\nread=%s\nbest_m=%i\nbest_i=%i\nbest_len=%i\n", READ_ID, READ, best_m, best_i, best_len) ;
				//exit(-1) ;
			}
		}
		
	}


	void find_poly(unsigned int &poly_length_start, unsigned int &poly_length_end, float frac=0.8)
	{
		int num_a_start = 0 ;
		int num_t_start = 0 ;
		int num_a_end = 0 ;
		int num_t_end = 0 ;
		poly_length_start = 0 ;
		poly_length_end = 0 ;
		
		for (unsigned int i=0; i<READ_LENGTH; i++)
		{
			if (READ[i]=='T' || READ[i]=='t')
				num_t_start++ ;
			if (READ[i]=='A' || READ[i]=='a')
				num_a_start++ ;
			if (READ[READ_LENGTH-i]=='A' || READ[READ_LENGTH-i]=='a')
				num_a_end++ ;
			if (READ[READ_LENGTH-i]=='T' || READ[READ_LENGTH-i]=='t')
				num_t_end++ ;
			if (((float)num_t_start)/i >= frac)
				poly_length_start = i ;
			if (((float)num_a_start)/i >= frac)
				poly_length_start = i ;
			if (((float)num_a_end)/i >= frac)
				poly_length_end = i ;
			if (((float)num_t_end)/i >= frac)
				poly_length_end = i ;
		}
	}

	bool is_full_poly(float frac=0.9)
		{
			int num_t=0, num_a=0 ;
			
			for (unsigned int i=0; i<READ_LENGTH; i++)
			{
				//TODO: dd sollte man nicht gleich ein toupper machen?
				if (READ[i]=='T' || READ[i]=='t')
					num_t++ ;
				if (READ[i]=='A' || READ[i]=='a')
					num_a++ ;
			}
			if (((float)num_a/READ_LENGTH)>=frac)
				return true ;
			if (((float)num_t/READ_LENGTH)>=frac)
				return true ;
			return false ;
		}
	
	Read const * get_orig() const { return orig_read ; } ;
//	void set_orig(Read* orig) { orig_read=orig ; } ;
	void printOn(FILE *file) const;

	int getNr() const {
		return _nr;
	}

	static int get_quality_offset() { return PRB_QUALITY_OFFSET ; } ; 
	static void set_quality_offset(int offset) { PRB_QUALITY_OFFSET=offset ; } ;
	
private:
	/**
	 * Copies data from another read - possibly truncated. This method is private,
	 * since the special case, src == this,  violates the const contract.
	 */
	void copyFrom(Read const &src, unsigned cutStart, unsigned cutEnd);

	int read_short_read();

	unsigned int READ_LENGTH;

	static int const _maxNrQualities = 3;
	char READ_QUALITY[_maxNrQualities][Config::MAX_READ_LENGTH + 1];
	char READ[Config::MAX_READ_LENGTH + 1];
	char READ_FORMAT;	// 0: fq, 1: fa, 2: flat
	char READ_ID[Config::MAX_READ_ID_LENGTH];
	int READ_PE_FLAG;
	int _nr;
	
	QueryFile &_queryFile;
	QueryFile::Location _location;

	Read const * orig_read ;

	static int PRB_QUALITY_OFFSET ;
};

inline std::ostream &operator<<(std::ostream &out, Read const &read) {
	read.printOn(out);
	return out;
}
