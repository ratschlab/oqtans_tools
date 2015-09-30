#pragma once

#include <palmapper/Genome.h>
#include <deque>
#include <string>
#include <stdlib.h> 

struct junction_str {
	int start; 
	int end; 
	int coverage;
	bool consensus ;
	char strand;
	std::string read_id ;
	std::string intron_string ;
	int junction_qual ;
};
typedef struct junction_str Junction;

struct exon_str {
	int start; 
	int end; 
	int chr;
	char strand;
};
typedef struct exon_str Exon;

class JunctionMap
{

public:
	JunctionMap(Genome const &genome_, int min_coverage_, 	std::vector<const char*> ACC_CONSENSUS_, std::vector<const char*> DON_CONSENSUS_, std::vector<const char*> ACC_CONSENSUS_REV_, std::vector<const char*> DON_CONSENSUS_REV_ ) ;
	~JunctionMap() ;

	void insert_junction(char strand, int chr, int start, int end, bool consensus, const char* intron_string, int junction_qual, const char* read_id, int coverage);
	int init_from_gffs(std::string &gff_fname);
	int report_to_gff(std::string &gff_fname);
	void filter_junctions();

	std::deque<Junction> * junctionlist ;

	void lock() 
	{ 
		pthread_mutex_lock( &junction_mutex) ; 
	}

	void unlock() 
	{ 
		pthread_mutex_unlock( &junction_mutex) ; 
	}
	

protected:
	bool is_consensus_intron(char strand, int chr, int start, int end);
	
	int init_from_gff(std::string &gff_fname);

	Genome const *genome;
	int min_coverage;

	pthread_mutex_t junction_mutex;
	std::vector<const char*> ACC_CONSENSUS, DON_CONSENSUS, ACC_CONSENSUS_REV, DON_CONSENSUS_REV ;	
   	


};

inline std::deque<Junction>::iterator  my_lower_bound ( std::deque<Junction>::iterator first, std::deque<Junction>::iterator  last, const int& value )
{
	std::deque<Junction>::iterator it;
	long int count, step;
	count = distance(first,last);

	while (count>0)
	{
		it = first;
		step = count/2;
		advance (it,step);
		
		if ( (*it).start < value) 
		{
			first=it; //++it;
			count-=step+1;
		}
		else count=step;
	}
	return first;
}
