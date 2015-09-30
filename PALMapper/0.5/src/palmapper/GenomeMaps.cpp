#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <zlib.h>
#include <limits>

#include "palmapper.h"
#include "dyn_prog/qpalma_dp.h"
#include <palmapper/Genome.h>
#include <palmapper/GenomeMaps.h>
#include <palmapper/Util.h>

clock_t GenomeMaps::last_report = 0 ;

GenomeMaps::GenomeMaps(Genome const &genome_)
{
	CHR_MAP_c = NULL ;
	CHR_MAP_i = NULL ;

	reported_repetitive_seeds = 0 ;
	reported_mapped_regions = 0 ;
	reported_mapped_reads = 0 ;
	reported_spliced_reads = 0 ;
    //reported_splice_sites = 0 ;

    covered_mapped_read_positions = 0 ;
    covered_mapped_read_positions_best = 0 ;
    covered_spliced_read_positions = 0 ;
    covered_spliced_read_positions_best = 0 ;
    covered_repetitive_seed_positions = 0 ;
    covered_repetitive_seed_positions_many1 = 0 ;
    covered_repetitive_seed_positions_many2 = 0 ;
    covered_mapped_region_positions = 0 ;
    //covered_splice_site_positions = 0 ;

	REPORT_REPETITIVE_SEED_DEPTH_EXTRA = 31 - MAX_INDEX_DEPTH ;

	genome = &genome_ ;

	if (_config.REPORT_REPETITIVE_SEEDS || _config.REPORT_MAPPED_REGIONS || _config.REPORT_MAPPED_READS || _config.REPORT_FILE!=NULL || _config.FILTER_BY_SPLICE_SITES || _config.QPALMA_USE_SPLICE_SITES)
	{
		init_reporting() ;

		if (!_config.REPORT_RESET)
		{
			read_reporting() ;
			do_reporting(1) ;
		}
	}

	if (_config.REPORT_GFF_FILE_NAME.size()>0)
		init_with_gff(_config.REPORT_GFF_FILE_NAME) ;
}

GenomeMaps::~GenomeMaps()
{
	clean_reporting() ;
}

#ifdef CHR_MAP_DNAARRAY
void GenomeMaps::to_dnaarray(int chr=-1)
{
	assert(CHR_MAP_c!=NULL) ;
OB	//fprintf(stdout, "to_dnaarray(%i)\n", chr) ;
	
	size_t start = 0 ;
	size_t end = 0 ;
	
	if (chr==-1)
	{
		start=0 ;
		end = genome->nrChromosomes() ;
	} 
	else
	{
		start=chr ;
		end=chr+1 ;
	}

	for (size_t i=start; i<end; i++) 
	{
		assert(CHR_MAP_c[i]!=NULL) ;

		for (size_t p=0; p<CHR_LENGTH[i]; p++)
		{
			unsigned char c = CHR_MAP_c[i][p] ;
#ifdef CHR_MAP_DNAARRAY
#ifndef CHR_MAP_DNAARRAY_2BIT
			c = c & (MASK_MAPPED_READ_BEST|MASK_MAPPED_READ|MASK_SPLICED_READ_BEST|MASK_SPLICED_READ) ;
#else
			// move the bit to the right position
			if (c & MASK_SPLICED_READ_BEST_ORIG)
				c = (c & MASK_MAPPED_READ_BEST) | MASK_SPLICED_READ_BEST ; 
			else
				c = c & MASK_MAPPED_READ_BEST ; 
#endif
#endif
			CHR_MAP_a[i]->set_elem(p, c) ;
		}

		delete[] CHR_MAP_c[i] ;
		CHR_MAP_c[i]=NULL ;
	}
}

void GenomeMaps::from_dnaarray(int chr = -1)
{
	assert(CHR_MAP_c!=NULL) ;
	//fprintf(stdout, "from_dnaarray(%i)\n", chr) ;

	size_t start = 0 ;
	size_t end = 0 ;
	
	if (chr==-1)
	{
		start=0 ;
		end = genome->nrChromosomes() ;
	} 
	else
	{
		start=chr ;
		end=chr+1 ;
	}

	for (size_t i=start; i<end; i++) 
	{
		assert(CHR_MAP_c[i]==NULL) ;
		
		try 
		{
			CHR_MAP_c[i] = new unsigned char[CHR_LENGTH[i]+1] ;
		}
		catch (std::bad_alloc&) 
		{
			fprintf(stderr, "[from_dnaarray]: ERROR: not enough memory for genome maps\n");
			exit(1);
		}
		
		for (size_t p=0; p<CHR_LENGTH[i]; p++)
		{
			CHR_MAP_c[i][p] = CHR_MAP_a[i]->get_elem(p) ;
#ifdef CHR_MAP_DNAARRAY
#ifndef CHR_MAP_DNAARRAY_2BIT
			assert((CHR_MAP_c[i][p] & (~(MASK_SPLICED_READ_BEST|MASK_MAPPED_READ_BEST|MASK_MAPPED_READ|MASK_SPLICED_READ)))==0) ;
#else
			assert((CHR_MAP_c[i][p] & (~(MASK_SPLICED_READ_BEST|MASK_MAPPED_READ_BEST)))==0) ;
			// move the bit to the right position
			if (CHR_MAP_c[i][p] & MASK_SPLICED_READ_BEST)
				CHR_MAP_c[i][p] = (CHR_MAP_c[i][p] & MASK_MAPPED_READ_BEST) | MASK_SPLICED_READ_BEST_ORIG ;
#endif
#endif
		}
		
	}
}
#endif

int GenomeMaps::init_reporting()
{
	try
	{
		CHR_MAP_c = new unsigned char*[genome->nrChromosomes()] ;
		if (_config.REPORT_GENOME_COVERAGE>0)
			CHR_MAP_i = new unsigned int*[genome->nrChromosomes()] ;
			
	}
	catch (std::bad_alloc&) 
	{
		fprintf(stderr, "ERROR : not enough memory for chr_map\n");
		exit(1);
	}

#ifdef CHR_MAP_DNAARRAY
	fprintf(stdout, "using compact map representation (class %s)\n", CHR_MAP_DNAARRAY_CLASS::get_class_name()) ;

	fprintf(stdout, "MASK_MAPPED_READ_BEST=%i\n", MASK_MAPPED_READ_BEST) ;
	fprintf(stdout, "MASK_MAPPED_READ=%i\n", MASK_MAPPED_READ) ;
	fprintf(stdout, "MASK_SPLICED_READ_BEST=%i\n", MASK_SPLICED_READ_BEST) ;
	fprintf(stdout, "MASK_SPLICED_READ=%i\n", MASK_SPLICED_READ) ;
#endif

	for (size_t i=0; i!=genome->nrChromosomes(); ++i) 
	{
		Chromosome const &chr = genome->chromosome(i);
#ifdef CHR_MAP_DNAARRAY
		CHR_MAP_a.push_back(new CHR_MAP_DNAARRAY_CLASS(CHR_LENGTH[i])) ;
		CHR_MAP_c[i]=NULL ;
		CHR_MAP_i[i]=NULL ;

		XXX
#else
		try 
		{
			CHR_MAP_c[i] = new unsigned char[chr.length()+1] ;
			if (_config.REPORT_GENOME_COVERAGE>0)
				CHR_MAP_i[i] = new unsigned int[chr.length()+1] ;
		}
		catch (std::bad_alloc&) 
		{
			fprintf(stderr, "ERROR : not enough memory for genome maps\n");
			exit(1);
		}
		memset(CHR_MAP_c[i], 0, chr.length()*sizeof(unsigned char)) ;
		if (_config.REPORT_GENOME_COVERAGE>0)
			memset(CHR_MAP_i[i], 0, chr.length()*sizeof(unsigned int)) ;
#endif
	}
	 
	covered_repetitive_seed_positions = 0 ;
	covered_repetitive_seed_positions_many1 = 0 ;
	covered_repetitive_seed_positions_many2 = 0 ;
	covered_mapped_region_positions = 0 ;
	covered_mapped_read_positions = 0 ;
	covered_mapped_read_positions_best = 0 ;
	covered_spliced_read_positions = 0 ;
	covered_spliced_read_positions_best = 0 ;
	//covered_splice_site_positions = 0 ;
	
	return 0 ;
}


int GenomeMaps::report_repetitive_seed(Chromosome const &chr, int chr_start, int count)
{
	//fprintf(stdout, "repetitive_seed:\tchr=%s\tpos=%i\tcount=%i\n", CHR_DESC[chr], chr_start, count) ;

	if (count<Config::REPORT_REPETITIVE_SEED_COUNT)
		return 0 ;
	if (!(chr_start>=0 && (size_t)chr_start<=chr.length()))
	{
		if (_config.VERBOSE>0)
			fprintf(stdout, "report_mapped_region: start out of bounds\n") ;
		if (chr_start<0)
			chr_start = 0 ;
		else
			chr_start = chr.length() ;
	}
	assert(chr_start>=0 && (size_t)chr_start<chr.length()) ;
	
	//fprintf(stdout, "repetitive_seed:\tchr=%s\tpos=%i\tcount=%i\n", CHR_DESC[chr], chr_start, count) ;
	reported_repetitive_seeds++ ;

	for (unsigned int i=chr_start; i<chr_start+_config.INDEX_DEPTH+REPORT_REPETITIVE_SEED_DEPTH_EXTRA && i<chr.length(); i++)
		if ((CHR_MAP(chr,i) & MASK_REPETITIVE_SEED) == 0)
		{
			covered_repetitive_seed_positions++ ;
			CHR_MAP_set(chr,i, CHR_MAP(chr,i) + MASK_REPETITIVE_SEED) ;
		}

	do_reporting() ;

	if (count<Config::REPORT_REPETITIVE_SEED_COUNT_MANY1)
		return 0 ;

	for (unsigned int i=chr_start; i<chr_start+_config.INDEX_DEPTH+REPORT_REPETITIVE_SEED_DEPTH_EXTRA && i<chr.length(); i++)
		if ((CHR_MAP(chr,i) & MASK_REPETITIVE_SEED_MANY1) == 0)
		{
			covered_repetitive_seed_positions_many1++ ;
			CHR_MAP_set(chr, i, CHR_MAP(chr,i) + MASK_REPETITIVE_SEED_MANY1) ;
		}

	if (count<Config::REPORT_REPETITIVE_SEED_COUNT_MANY2)
		return 0 ;

	for (unsigned int i=chr_start; i<chr_start+_config.INDEX_DEPTH+REPORT_REPETITIVE_SEED_DEPTH_EXTRA && i<chr.length(); i++)
		if ((CHR_MAP(chr,i) & MASK_REPETITIVE_SEED_MANY2) == 0)
		{
			covered_repetitive_seed_positions_many2++ ;
			CHR_MAP_set(chr, i, CHR_MAP(chr,i) + MASK_REPETITIVE_SEED_MANY2) ;
		}


	return 1 ;
}

int GenomeMaps::report_mapped_region(Chromosome const &chr, int chr_start, int chr_end, int num_matches)
{
	reported_mapped_regions++ ;
	//fprintf(stdout, "mapped_region:\tchr=%s\tstart=%i\tend=%i\tmatches=%i\n", CHR_DESC[chr], chr_start, chr_end, num_matches) ;
	//assert(chr_start>=0 && chr_start<CHR_LENGTH[chr]) ;
	//assert(chr_end>=0 && chr_end<CHR_LENGTH[chr]) ;

	if (!(chr_start>=0 && (size_t)chr_start<=chr.length()))
	{
		if (_config.VERBOSE>0)
			fprintf(stdout, "report_mapped_region: start out of bounds\n") ;
		if (chr_start<0)
			chr_start = 0 ;
		else
			chr_start = chr.length() ;
	}
	if (!(chr_end>=0 && (size_t)chr_end<=chr.length()))
	{
		if (_config.VERBOSE>0)
			fprintf(stdout, "report_mapped_region: end out of bounds\n") ;
		if (chr_end<0)
			chr_end = 0 ;
		else
			chr_end = chr.length() ;
	}

	for (int i=chr_start; i<chr_end; i++)
		if ((CHR_MAP(chr,i) & MASK_MAPPED_REGION) ==0 )
		{
			covered_mapped_region_positions++ ;
			CHR_MAP_set(chr, i, CHR_MAP(chr,i) + MASK_MAPPED_REGION) ;
		}
	
	return 0 ;
}

int GenomeMaps::report_mapped_read(Chromosome const &chr, int start, int end, int num_matches, int nbest_hit)
{
	reported_mapped_reads++ ;
	
	if (!(start>=0 && (size_t)start<=chr.length()))
	{
		if (_config.VERBOSE>0)
			fprintf(stdout, "report_mapped_read: start out of bounds (start=%i no in [0,%i))\n", start, chr.length()) ;
		if (start<0)
			start = 0 ;
		else
			start = chr.length() ;
	}
	
	if (!(end>=0 && (size_t)end<=chr.length()))
	{
		if (_config.VERBOSE>0)
			fprintf(stdout, "report_mapped_read: end out of bounds (end=%i not in [0,%i))\n", end, chr.length()) ;
		if (end<0)
			end = 0 ;
		else
			end = chr.length() ;
	}

	//fprintf(stdout, "mapped_read:  \tchr=%s\tstart=%i\tend=%i\tmatches=%i\tnbest=%i\n", chr.desc(), start, end, num_matches, nbest_hit) ;
	int start_map=start-Config::QPALMA_USE_MAP_WINDOW ;
	if (start_map<0)
		start_map=0 ;
	
	int end_map=end+Config::QPALMA_USE_MAP_WINDOW ;
	if (end_map>(int)chr.length())
		end_map=chr.length() ;
	
	for (int i=start_map; i<end_map; i++)
	{
		if (nbest_hit==0 && (CHR_MAP(chr,i) & MASK_MAPPED_READ_BEST) == 0 )
		{
			covered_mapped_read_positions_best++ ;
			CHR_MAP_set(chr, i, CHR_MAP(chr,i)+MASK_MAPPED_READ_BEST) ;
		}
#ifndef CHR_MAP_DNAARRAY
		if ((CHR_MAP(chr,i) & MASK_MAPPED_READ) == 0 )
		{
			covered_mapped_read_positions++ ;
			CHR_MAP_set(chr, i, CHR_MAP(chr,i)+ MASK_MAPPED_READ) ;
		}
#else
#ifndef CHR_MAP_DNAARRAY_2BIT
		if ((CHR_MAP(chr,i) & MASK_MAPPED_READ) == 0 )
		{
			covered_mapped_read_positions++ ;
			CHR_MAP_set(chr, i, CHR_MAP(chr,i)+ MASK_MAPPED_READ) ;
		}
#endif
#endif
	}
	
	if (_config.REPORT_GENOME_COVERAGE>0)
		for (int i=start; i<end; i++)
		{
			if (CHR_MAP_cov(chr,i)<std::numeric_limits<unsigned int>::max())
			{
				CHR_MAP_set_cov(chr, i, CHR_MAP_cov(chr,i)+1) ;
			}
		}
	
	do_reporting() ;
				
	return 0 ;
}

int GenomeMaps::report_spliced_read(Chromosome const &chr, std::vector<int> & exons, int num_matches, int nbest_hit)
{
	reported_spliced_reads++ ;

	//fprintf(stdout, "mapped_read:  \tchr=%s\tstart=%i\tend=%i\tmatches=%i\tnbest=%i\n", CHR_DESC[chr], start, end, num_matches, nbest_hit) ;
	
	if (exons.size()==2)
	{
		return report_mapped_read(chr, exons[0], exons[1], num_matches, nbest_hit) ;
	} 

	for (size_t e=0; e<exons.size(); e+=2)
	{
		int start = exons[e] ;
		int end = exons[e+1] ;
		
		if (!(start>=0 && (size_t)start<=chr.length()))
		{
			if (_config.VERBOSE>0)
				fprintf(stdout, "report_spliced_read: start out of bounds\n") ;
			if (start<0)
				start = 0 ;
			else
				start = chr.length();
		}
		
		if (!(end>=0 && (size_t)end<=chr.length()))
		{
			if (_config.VERBOSE>0)
				fprintf(stdout, "report_spliced_read: end out of bounds\n") ;
			if (end<0)
				end = 0 ;
			else
				end = chr.length();
		}
		
		//assert(start>=0 && (size_t)start<=CHR_LENGTH[chr]) ;
		//assert(end>=0 && (size_t)end<=CHR_LENGTH[chr]) ;
		int start_map = start - Config::QPALMA_USE_MAP_WINDOW ;
		if (start_map<0)
			start_map=0 ;
		int end_map = end + Config::QPALMA_USE_MAP_WINDOW ;
		if (end_map>(int)chr.length())
			end_map=chr.length() ;
		
		for (int i=start_map; i<end_map; i++)
		{
			if (nbest_hit==0 && (CHR_MAP(chr,i) & MASK_SPLICED_READ_BEST) == 0 )
			{
				covered_spliced_read_positions_best++ ;
				CHR_MAP_set(chr, i, CHR_MAP(chr,i)+MASK_SPLICED_READ_BEST) ;
			}
#ifndef CHR_MAP_DNAARRAY
			if ((CHR_MAP(chr,i) & MASK_SPLICED_READ) == 0 )
			{
				covered_spliced_read_positions++ ;
				CHR_MAP_set(chr, i, CHR_MAP(chr,i)+MASK_SPLICED_READ) ;
			}
#else
#ifndef CHR_MAP_DNAARRAY_2BIT
			if ((CHR_MAP(chr,i) & MASK_SPLICED_READ) == 0 )
			{
				covered_spliced_read_positions++ ;
				CHR_MAP_set(chr, i, CHR_MAP(chr,i)+MASK_SPLICED_READ) ;
			}
#else 
			assert(MASK_SPLICED_READ_BEST<4) ;
#endif
#endif
		}
		if (_config.REPORT_GENOME_COVERAGE>0)
			for (int i=start; i<end; i++)
			{
				if (CHR_MAP_cov(chr,i)<std::numeric_limits<unsigned int>::max())
				{
					CHR_MAP_set_cov(chr, i, CHR_MAP_cov(chr,i)+1) ;
				}
			}
	}
						
	return 0 ;
}

/*
int GenomeMaps::report_splice_site(int chr, int pos, char strand, char type)
{
	reported_splice_sites++ ;

	assert(pos>=0 && pos<CHR_LENGTH[chr]) ;

	int MASK_SPLICE_SITE = (strand == '+') ? MASK_SPLICE_SITE_P : MASK_SPLICE_SITE_N ;
	
	if ( (CHR_MAP[chr][pos] & MASK_SPLICE_SITE) == 0 )
	{
		covered_splice_site_positions++ ;
		CHR_MAP[chr][pos] += MASK_SPLICE_SITE ;
	}
	// this storage scheme assumes that acceptor and donor splice sites do not occur at the same position
    // (I'm running out of bits here ... )

	// acc = true
	if ((type == 'a') && ((CHR_MAP[chr][pos] & MASK_SPLICE_SITE_ACC) == 0) )
	{
		CHR_MAP[chr][pos] += MASK_SPLICE_SITE_ACC ;
	}
	// don = false
	if ((type == 'd') && ((CHR_MAP[chr][pos] & MASK_SPLICE_SITE_ACC) != 0) )
	{
		CHR_MAP[chr][pos] -= MASK_SPLICE_SITE_ACC ;
	}
	assert(type=='a' || type=='d') ;
	
	return 0 ;
}
*/

int GenomeMaps::do_reporting(int force)
{
	if (force || (clock()-last_report)/CLOCKS_PER_SEC>10)
	{
		last_report=clock() ;

		fprintf(stdout, "\n") ;
		if (_config.VERBOSE>0)
			fprintf(stdout, "[report] repetitive_seeds (%i)\t=\t %8i\t(%8i, %8i, %8i positions)\n", REPORT_REPETITIVE_SEED_DEPTH_EXTRA+_config.INDEX_DEPTH, reported_repetitive_seeds, covered_repetitive_seed_positions, covered_repetitive_seed_positions_many1, covered_repetitive_seed_positions_many2) ;
		if (_config.VERBOSE>0)
			fprintf(stdout, "[report] mapped_regions   \t=\t %8i\t(%8i positions)\n", reported_mapped_regions, covered_mapped_region_positions) ;
		if (_config.VERBOSE>0)
			fprintf(stdout, "[report] mapped_reads     \t=\t %8i\t(%8i, %8i positions)\n", reported_mapped_reads, covered_mapped_read_positions_best, covered_mapped_read_positions) ;
		if (_config.VERBOSE>0)
			fprintf(stdout, "[report] spliced_reads     \t=\t %8i\t(%8i, %8i positions)\n", reported_spliced_reads, covered_spliced_read_positions_best, covered_spliced_read_positions) ;
		//fprintf(stderr, "[report] splice_sites     =\t %8i\t(%8i positions)\n\n", reported_splice_sites, covered_splice_site_positions) ;
	}

	return 0 ;
}

int GenomeMaps::read_reporting()
{
	const char* fname = _config.REPORT_FILE ;
	if (fname==NULL || strlen(_config.REPORT_FILE)==0)
		return -1 ;

	gzFile fd = gzopen(fname, "rb") ;
	if (!fd)
		return -1 ;

	covered_repetitive_seed_positions = 0 ;
	covered_repetitive_seed_positions_many1 = 0 ;
	covered_repetitive_seed_positions_many2 = 0 ;
	covered_mapped_region_positions = 0 ;
	covered_mapped_read_positions = 0 ;
	covered_mapped_read_positions_best = 0 ;
	covered_spliced_read_positions = 0 ;
	covered_spliced_read_positions_best = 0 ;
	//covered_splice_site_positions = 0 ;

    // sanity checking
#ifdef CHR_MAP_DNAARRAY
#ifndef CHR_MAP_DNAARRAY_2BIT
	assert(strcmp(CDNAArray::get_class_name(), CHR_MAP_DNAARRAY_CLASS::get_class_name())==0) ;
#else
	assert(strcmp(CDNAArray4::get_class_name(), CHR_MAP_DNAARRAY_CLASS::get_class_name())==0) ;
#endif
#endif

	for (size_t i=0; i<genome->nrChromosomes(); i++)
	{
		Chromosome const &chr = genome->chromosome(i);
#ifdef CHR_MAP_DNAARRAY
		bool to_be_transfered=false ;
		if (CHR_MAP_c[i]==NULL)
		{
			from_dnaarray(i) ;
			to_be_transfered = true ;
		}
#endif
		assert(CHR_MAP_c[i]!=NULL) ;
		
		int ret=gzread(fd, CHR_MAP_c[i], chr.length()) ;
		if ((size_t)ret!=chr.length())
		{
			fprintf(stderr, "error: Genome map inconsistent (%s)\n", fname) ;
			exit(-2) ;
		}

#ifdef CHR_MAP_DNAARRAY
		if (to_be_transfered)
			to_dnaarray(i) ;
#endif
		
		for (size_t j=0; j<chr.length(); j++)
		{
            // bits are reordered for CDNAArray4
			if ((CHR_MAP(chr,j) & MASK_MAPPED_READ_BEST)!=0)
				covered_mapped_read_positions_best++ ;
			if ((CHR_MAP(chr,j) & MASK_MAPPED_READ)!=0)
				covered_mapped_read_positions++ ;
			if ((CHR_MAP(chr,j) & MASK_SPLICED_READ_BEST)!=0)
				covered_spliced_read_positions_best++ ;
			if ((CHR_MAP(chr,j) & MASK_SPLICED_READ)!=0)
				covered_spliced_read_positions++ ;
			if ((CHR_MAP(chr,j) & MASK_MAPPED_REGION)!=0)
				covered_mapped_region_positions++ ;
			if ((CHR_MAP(chr,j) & MASK_REPETITIVE_SEED)!=0)
				covered_repetitive_seed_positions++ ;
			if ((CHR_MAP(chr,j) & MASK_REPETITIVE_SEED_MANY1)!=0)
				covered_repetitive_seed_positions_many1++ ;
			if ((CHR_MAP(chr,j) & MASK_REPETITIVE_SEED_MANY2)!=0)
				covered_repetitive_seed_positions_many2++ ;
		}
	}
	gzclose(fd) ;

	fprintf(stdout, "read QPALMA reporting map from %s\n", fname) ;
	return 0 ;
}

int GenomeMaps::write_reporting()
{
	const char* fname = _config.REPORT_FILE ;
	if (fname==NULL || strlen(_config.REPORT_FILE)==0)
		return -1 ;
	if (_config.REPORT_FILE_READONLY)
		return -1 ;
	
	gzFile fd = gzopen(fname, "wb6") ;
	if (!fd)
		return -1 ;
	fprintf(stdout, "writing QPALMA reporting map to %s:\n", fname) ;
	for (size_t i=0; i<genome->nrChromosomes(); i++)
	{
		Chromosome const &chr = genome->chromosome(i);
#ifdef CHR_MAP_DNAARRAY
		bool to_be_deleted=false ;
		if (CHR_MAP_c[i]==NULL)
		{
			from_dnaarray(i) ; // this may move bits around
			to_be_deleted=true ;
		}
#endif		
		int ret = gzwrite(fd, CHR_MAP_c[i], chr.length()) ;
		assert((size_t)ret==chr.length()) ;
		int num_covered = 0 ;
		for (size_t p=0; p<chr.length(); p++)
			if (CHR_MAP_c[i][p])
				num_covered++ ;
		if (_config.VERBOSE>0)
			fprintf(stdout, "  chr\t%s\thas \t%i\t/\t%i\t(%2.1f%%) covered positions\n", chr.desc(), num_covered, chr.length(), 100.0*num_covered/chr.length()) ;

#ifdef CHR_MAP_DNAARRAY
		if (to_be_deleted)
		{
			delete[] CHR_MAP_c[i] ;
			CHR_MAP_c[i]=NULL ;
		}
#endif
	}
	gzclose(fd) ;

	if (_config.VERBOSE>0)
		fprintf(stdout, "done.\n") ;

	
	return 0 ;
}

int GenomeMaps::write_cov_reporting()
{
	const char *fname = _config.REPORT_GENOME_COVERAGE_FILE ;
	
	if (fname==NULL || strlen(fname)==0)
		return -1 ;
	
	if (_config.REPORT_GENOME_COVERAGE==1) // binary
	{
		gzFile fd = gzopen(fname, "wb6") ;
		if (!fd)
			return -1 ;
		fprintf(stdout, "writing coverage map in binary format to %s:\n", fname) ;
		for (size_t i=0; i<genome->nrChromosomes(); i++)
		{
			Chromosome const &chr = genome->chromosome(i);
			int ret = gzwrite(fd, CHR_MAP_i[i], chr.length()*sizeof(unsigned int)) ;
			assert((size_t)ret==chr.length()*sizeof(unsigned int)) ;
			int num_covered = 0 ;
			for (size_t p=0; p<chr.length(); p++)
				if (CHR_MAP_i[i][p]>0)
					num_covered++ ;
			if (_config.VERBOSE>0)
				fprintf(stdout, "  chr\t%s\thas \t%i\t/\t%i\t(%2.1f%%) covered positions\n", chr.desc(), num_covered, chr.length(), 100.0*num_covered/chr.length()) ;
		}
		gzclose(fd) ;
	}

	if (_config.REPORT_GENOME_COVERAGE==2) // wiggle
	{
		FILE *fd = fopen(fname, "w+") ;
		if (!fd)
			return -1 ;
		fprintf(stdout, "writing coverage map in wiggle format to %s:\n", fname) ;
		fprintf(fd, "track type=wiggle_0 name=read_coverage description=read_coverage\n") ;
		for (size_t i=0; i<genome->nrChromosomes(); i++)
		{
			Chromosome const &chr = genome->chromosome(i);
			fprintf(fd, "variableStep chrom=%s\n", chr.desc()) ;
			int num_covered = 0 ;
			for (unsigned int j=0; j<chr.length(); j++)
			{
				if (CHR_MAP_i[i][j]>0)
				{
					fprintf(fd, "%i %i\n", j, CHR_MAP_i[i][j]) ;
					num_covered++ ;
				}
			}
			if (_config.VERBOSE>0)
				fprintf(stdout, "  chr\t%s\thas \t%i\t/\t%i\t(%2.1f%%) covered positions\n", chr.desc(), num_covered, chr.length(), 100.0*num_covered/chr.length()) ;
		}
		fclose(fd) ;
	}

	if (_config.VERBOSE>0)
		fprintf(stdout, "done.\n") ;
	return 0 ;
}

int GenomeMaps::clean_reporting()
{
	if (CHR_MAP_c!=NULL)
	{
		for (size_t i=0; i!=genome->nrChromosomes(); ++i) 
			delete[] CHR_MAP_c[i] ;
		
		delete[] CHR_MAP_c ;
		CHR_MAP_c=NULL ;
	}
	if (CHR_MAP_i!=NULL)
	{
		for (size_t i=0; i!=genome->nrChromosomes(); ++i) 
			delete[] CHR_MAP_i[i] ;
		
		delete[] CHR_MAP_i ;
		CHR_MAP_i=NULL ;
	}

#ifdef CHR_MAP_DNAARRAY
	for (size_t i=0; i<CHR_MAP_a.size(); i++)
		delete CHR_MAP_a[i] ;
#endif

	return 0 ;
}

int GenomeMaps::init_with_gff(std::string &gff_fname)
{
	fprintf(stdout, "initializing genome map with GFF file %s\n", gff_fname.c_str()) ;

	FILE * fd=Util::openFile(gff_fname.c_str(), "r") ;
	if (!fd)
		return -1 ;
	int exon_lines=0 ;
	
	while (!feof(fd))
	{
		char chr_name[1000], source[1000], type[1000], properties[1000], strand, tmp1, tmp2 ;
		int start, end ;

		Util::skip_comment_lines(fd) ;
		
		int num = fscanf(fd, "%1000s\t%1000s\t%1000s\t%i\t%i\t%c\t%c\t%c\t%1000s\n", chr_name, source, type, &start, &end, &tmp1, &strand, &tmp2, properties) ;  
		if (num!=9)
		{
			if (feof(fd))
				break ;
			fprintf(stdout, "gff line only contained %i columns, aborting\n", num) ;
		}
		
		if (strcmp(type, "exon")==0)
		{
			exon_lines++ ;
			
			int chr_idx = genome->find_desc(chr_name) ;
			if (chr_idx==-1)
			{
				fprintf(stderr, "chromosome %s not found. known chromosome names:\n", chr_name) ;
				genome->print_desc(stderr) ;
				return -1 ;
			}
			Chromosome const &chr = genome->chromosome(chr_idx) ;
			
			if (!(start>=0 && (size_t)start<=chr.length()))
			{
				fprintf(stderr, "init_with_gff: start out of bounds\n") ;
				if (start<0)
					start = 0 ;
				else
					start = chr.length() ;
			}
			
			if (!(end>=0 && (size_t)end<=chr.length()))
			{
				fprintf(stderr, "init_with_gff: end out of bounds\n") ;
				if (end<0)
					end = 0 ;
				else
					end = chr.length() ;
			}

			for (int i=start; i<end; i++)
			{
				if ( (CHR_MAP(chr,i) & MASK_MAPPED_READ_BEST) == 0 )
				{
					CHR_MAP_set(chr, i, CHR_MAP(chr,i)+MASK_MAPPED_READ_BEST) ;
				}
#ifndef CHR_MAP_DNAARRAY
				if ((CHR_MAP(chr,i) & MASK_MAPPED_READ) == 0 )
				{
					CHR_MAP_set(chr, i, CHR_MAP(chr,i)+ MASK_MAPPED_READ) ;
				}
#else
#ifndef CHR_MAP_DNAARRAY_2BIT
				if ((CHR_MAP(chr,i) & MASK_MAPPED_READ) == 0 )
				{
					CHR_MAP_set(chr, i, CHR_MAP(chr,i)+ MASK_MAPPED_READ) ;
				}
#endif
#endif
				
			}
		}
	}
	fclose(fd) ;

	fprintf(stdout, "read %i exon lines\n", exon_lines) ;
	
	return 0 ;
}
