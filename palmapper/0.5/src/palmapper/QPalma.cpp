// Authors: Gunnar Raetsch, Geraldine Jean, Lisa Thalheim
// Copyright (C) 2009-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include <assert.h>
#include <string>
#include <sys/time.h>
#include <time.h>

#include "palmapper.h"
#include "IntervalQuery.h"
#include "dyn_prog/qpalma_dp.h"

#include <palmapper/Genome.h>
#include <palmapper/QPalma.h>
#include <palmapper/JunctionMap.h>

#define MAX_NUM_LONG_HITS _config.SPLICED_MAX_NUM_ALIGNMENTS 
#define MAX_MAP_REGION_SIZE _config.QPALMA_USE_MAP_MAX_SIZE


void get_annotated_splice_positions( std::vector<int> &pos, JunctionMap &annotatedjunctions,const char * type, int start, int end, int chr, char strand)
{
	// find a lower bound on the index with binary search
	annotatedjunctions.lock() ;
	//std::deque<Junction>::iterator it = annotatedjunctions.junctionlist[chrN].begin();
	std::deque<Junction>::iterator it = my_lower_bound(annotatedjunctions.junctionlist[chr].begin(), annotatedjunctions.junctionlist[chr].end(), start-100) ;

	while (it != annotatedjunctions.junctionlist[chr].end() && (*it).start <= end)
	{

		if ((*it).strand!=strand){
			it++;
			continue;
		}
		//fprintf(stdout,"Splice junction: %i-%i %c (%i-%i %c)\n", (*it).start, (*it).end, (*it).strand,start, end,strand );		
		if ((*it).start>=start && (*it).start< end){
			if (strcmp(type, "don")==0 && strand=='+') 
				pos.push_back((*it).start);
			if (strcmp(type, "acc")==0 && strand=='-') 
				pos.push_back((*it).start);
		}
		if ((*it).end>=start && (*it).end< end){
			if (strcmp(type, "acc")==0 && strand=='+') 
				pos.push_back((*it).end);
			if (strcmp(type, "don")==0 && strand=='-') 
				pos.push_back((*it).end);
		}
		it++;
		
	}
	
	annotatedjunctions.unlock();
	
	
}

void get_annotated_splice_sites(std::vector<int> &acc_pos, std::vector<int> &don_pos, JunctionMap &annotatedjunctions, int start, int end, int chr, char strand)
{
	//fprintf(stdout,"Get splice sites\n");
	// find a lower bound on the index with binary search
	annotatedjunctions.lock() ;
	//std::deque<Junction>::iterator it = annotatedjunctions.junctionlist[chrN].begin();
	std::deque<Junction>::iterator it = my_lower_bound(annotatedjunctions.junctionlist[chr].begin(), annotatedjunctions.junctionlist[chr].end(), start-100) ;

	while (it != annotatedjunctions.junctionlist[chr].end() && (*it).start <= end)
	{


		if ((*it).strand!=strand){
			it++;
			continue;
		}

		//fprintf(stdout,"Splice junction: %i-%i %c (%i-%i %c) cons=%i\n", (*it).start, (*it).end, (*it).strand,start, end,strand, (*it).consensus );		
		
		if ((*it).start>=start && (*it).start< end){
			if (strand=='+'){
				don_pos.push_back((*it).start);
			}
			else
				acc_pos.push_back((*it).start);
		}
		if ((*it).end>=start && (*it).end< end){
			if (strand=='+'){
				acc_pos.push_back((*it).end);
			}
			
			else
				don_pos.push_back((*it).end);
		}
		it++;
		

			
	}
	
	annotatedjunctions.unlock();
}


int QPalma::get_transcription_direction(int side,int orientation) const
{

	//Return 1 if the transcription is on the forward direction
	//Return -1 if the transcription is on the negative direction
	//Return 0 if we don't know

	//No strand specific information: need to try all combinations to look for splice sites
	if (side<=-1 || _config.PROTOCOL == PROTOCOL_UNSTRANDED)
		return 0;
	
	
	//Left reads: 
	//same orientation than the transcription for first strand protocol
	//opposite orientation than the transcription for second strand protocol
	if (side==1){
		if (_config.PROTOCOL == PROTOCOL_FIRST)
			return (orientation == 0)?1:-1;
		else
			return (orientation == 0)?-1:1;
	}
	
	//Right reads:
	//opposite orientation than the transcription for first strand protocol
	//same orientation than the transcription for second strand protocol
	if (_config.PROTOCOL == PROTOCOL_FIRST)
		return (orientation == 0)?-1:1;
	else
		return (orientation == 0)?1:-1;
		   
}

QPalma::QPalma(Genome* genome_, GenomeMaps* genomemaps_, int verbosity_)
//  : 	verbosity(3), MIN_NUM_MATCHES(_config.QPALMA_MIN_NUM_MATCHES)
	: 	verbosity(verbosity_), MIN_NUM_MATCHES(_config.QPALMA_MIN_NUM_MATCHES)
{
	genome=genome_ ;
	genomemaps = genomemaps_ ;
	alignment_parameters = NULL;
	
	if (_config.SPLICED_HITS)
	{
		if (_config.QPALMA_FILE.length()>0 || _config.NO_QPALMA)
		{
			int ret=init_alignment_parameters(_config.QPALMA_FILE) ;
			if (ret!=0)
			{
				fprintf(stderr, "init_alignment_parameters failed\n") ;
				exit(ret) ;
			}
		}
		else
		{
				fprintf(stderr, "alignment parameter file not provided (use -qpalma option)\n") ;
				exit(-1) ;
		}
		
		if (_config.ACC_FILES.length()>0)
		{
			int ret=check_splice_files(_config.ACC_FILES) ;
			if (ret!=0)
			{
				fprintf(stderr, "check_splice_file failed\n") ;
				exit(ret) ;
			}
		}
		else
		{
			if (!_config.NO_SPLICE_PREDICTIONS && !_config.ANNOTATED_SPLICE_SITES_FILE.length()>0)
			{
				fprintf(stderr, "acceptor splice site predictions file not provided (use -acc option)\n") ;
				exit(-1) ;
			}
		}
		
		if (_config.DON_FILES.length()>0)
		{
			int ret=check_splice_files(_config.DON_FILES) ;
			if (ret!=0)
			{
				fprintf(stderr, "check_splice_file failed\n") ;
				exit(ret) ;
			}
		}
		else
		{
			if (!_config.NO_SPLICE_PREDICTIONS && !_config.ANNOTATED_SPLICE_SITES_FILE.length()>0)
			{
				fprintf(stderr, "donor splice site predictions file not provided (use -acc option)\n") ;
				exit(-1) ;
			}
		}
	}
}

QPalma::Result::Result(Read const &read, QPalma const &qpalma)
	:	_qpalma(qpalma), _read(read)
{
	if (&qpalma == NULL)
		return;
	qpalma_filter_reason = -1 ;
	for (int ori=0; ori<2; ori++)
		for (uint32_t i = 0; i < qpalma.genome->nrChromosomes(); i++)
		{
			std::vector<region_t *> r; 
			regions[ori].push_back(r);
		}
}


QPalma::~QPalma()
{
	clean_alignment_parameters() ;
}


int QPalma::check_splice_files(std::string file_template)
{
	char basename[1000] ;
	
	for (int chr=1; chr<(int)genome->nrChromosomes(); chr++)
	{
		char strand='+' ;
		
		while (1)
		{
			sprintf(basename, file_template.c_str(), chr, strand) ;
		    char posname[1000] ;
			sprintf(posname, "%s.pos", basename) ;
		    char confcumname[1000] ;
			sprintf(confcumname, "%s.Conf_cum", basename) ;

			FILE *fd=fopen(posname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", posname) ;
				return -1 ;
			}
			int ret = fseek(fd, -4, SEEK_END) ;
			if (ret==0) // file may be empty
			{
				int buf;
				int retu = fread(&buf, sizeof(int), 1, fd) ;
				if (retu!=1)
				{
					fprintf(stderr, "fread error on %s\n", posname) ;
					return -1 ;
				}
				if (buf>(int)genome->chromosome(chr-1).length())
				{
					fprintf(stderr, "splice site position (%i) greater than chromosome length (%i)\n", buf, genome->chromosome(chr-1).length()) ;
					return -1 ;
				}
				//else
				//	fprintf(stdout, "last_pos=%i, size=%i\n", buf, genome->chromosome(chr-1).length()) ;
			}
			//else
			//fprintf(stderr, "fseek %s failed\n", posname) ;
			fclose(fd) ;

			fd=fopen(confcumname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", confcumname) ;
				return -1 ;
			}
			fclose(fd) ;
			
			struct stat pos_stat ;
			int statret = stat(posname, & pos_stat) ;
			if (statret!=0)
			{
			    fprintf(stderr, "stat on %s failed\n", posname) ;
			    return -1 ;
			}

			struct stat confcum_stat ;
			statret = stat(confcumname, & confcum_stat) ;
			if (statret!=0)
			{
			    fprintf(stderr, "stat on %s failed\n", confcumname) ;
			    return -1 ;
			}
			
			if (pos_stat.st_size!=confcum_stat.st_size)
			    fprintf(stdout, "**The two files %s and %s have different sizes (%ld, %ld).**\n**Continuing anyway.**\n\n", posname, confcumname, (long int)pos_stat.st_size, (long int)confcum_stat.st_size);

			if (strand=='+')
				strand='-' ;
			else
				break ;
		}
	}
	return 0 ;
}

int QPalma::compare_double(const void *a, const void *b) 
{
	double *ia = (double *) a;
	double *ib = (double *) b;

	return (int) (*ia - *ib);
}

int QPalma::map_splice_sites(std::string file_template, char type, float &splice_site_threshold, bool estimate_thresh, bool do_report)
{
	char basename[1000] ;
	
	for (int chr=0; chr<(int)genome->nrChromosomes() && (chr==0 || do_report); chr++)
	{
		char strand='+' ;
		
		while (1)
		{
			sprintf(basename, file_template.c_str(), chr+1, strand) ;
		    char posname[1000] ;
			sprintf(posname, "%s.pos", basename) ;
		    char confcumname[1000] ;
			sprintf(confcumname, "%s.Conf_cum", basename) ;

			FILE *fd=fopen(posname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", posname) ;
				return -1 ;
			}
			fclose(fd) ;

			fd=fopen(confcumname, "r") ;
			if (fd==NULL)
			{
				fprintf(stderr, "%s does not exist\n", confcumname) ;
				return -1 ;
			}
			fclose(fd) ;
			
			{
				IntervalQuery iq;
				const int num_scores = 1;
				const char *score_names[num_scores] = { "Conf_cum" };
				const int num_intervals = 1 ;
				int interval_matrix[num_intervals * 2];
				int cum_length[num_intervals + 1];
				cum_length[0] = 0;
				interval_matrix[0] = 0 ;
				interval_matrix[1] = genome->chromosome(chr).length();
				
				// acc
				sprintf(basename, file_template.c_str(), chr+1, strand) ;
				
				int acc_size = iq.query(basename, (char**) score_names, num_scores,
										interval_matrix, num_intervals);
				if (acc_size < 0) 
					return -1;
				
				int* acc_pos = NULL ;
				int* acc_index = NULL ;
				double* acc_score = NULL ;
				try
				{
					acc_pos = new int[acc_size] ;
					acc_index = new int[acc_size] ;
					acc_score = new double[acc_size] ;
				}
				catch (std::bad_alloc&)
				{
					fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
					delete[] acc_pos ;
					delete[] acc_index ;
					delete[] acc_score ;
					return -1;
				}
				
				iq.getResults(acc_pos, acc_index, acc_score);
				iq.cleanup();

				if (strand=='+' && chr==0 && estimate_thresh)
				{
					std::vector<double> acc_score2 ;
					
					try
					{
						for (int i=0; i<acc_size; i++)
							acc_score2.push_back(acc_score[i]) ;
					}
					catch (std::bad_alloc&)
					{
						fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
						delete[] acc_pos ;
						delete[] acc_index ;
						delete[] acc_score ;
						return -1;
					}
					
					std::sort(acc_score2.begin(), acc_score2.end()) ; 
					
					if (_config.VERBOSE)
						fprintf(stdout, "estimate_threshold: min=%1.3f  max=%1.3f %1.2f%% percentile=%1.3f\n", acc_score2[0], acc_score2[acc_size-1], splice_site_threshold*100, acc_score2[acc_size-1-(int)(acc_size*splice_site_threshold)]) ;
					
					splice_site_threshold=acc_score2[acc_size-1-(int)(acc_size*splice_site_threshold)] ;
				}

				/*int num_reported = 0 ;
				  if (do_report)
				  for (int i=0; i<acc_size; i++)
				  {
				  if (acc_score[i]>splice_site_threshold)
				  {
				  num_reported++ ;
				  report_splice_site(chr, acc_pos[i], strand, type) ;
				  }
				  }*/
				
				delete[] acc_pos ;
				delete[] acc_index ;
				delete[] acc_score ;
				
				//fprintf(stderr, "reported %i splice sites for %i %c\n", num_reported, chr, strand) ;
			}
			
			if (strand=='+')
				strand='-' ;
			else
				break ;
		}
	}

	return 0 ;
	
}

int QPalma::init_alignment_parameters(std::string qpalma_file)
{
	// initialize parameters for alignment
	if (alignment_parameters == NULL)
	{
		alignment_parameters = (struct alignment_parameter_struct*) malloc(
			sizeof(struct alignment_parameter_struct));
		
		if (alignment_parameters == NULL) 
		{
			fprintf(stderr, "[init_alignment_parameters] Could not allocate memory\n");
			return -1;
		}
		
		int ret = init_spliced_align(qpalma_file.c_str(),
									 alignment_parameters->h, alignment_parameters->a,
									 alignment_parameters->d, alignment_parameters->qualityPlifs,
									 alignment_parameters->num_qualityPlifs,
									 alignment_parameters->matchmatrix,
									 alignment_parameters->matchmatrix_dim,
									 alignment_parameters->quality_offset);
		if (ret < 0)
			return ret;
		return 0 ;
	} 
	else
		assert(0) ;
	return -1 ;
}

int QPalma::clean_alignment_parameters()
{
	if (alignment_parameters==NULL)
		return -1 ;
	
	free(alignment_parameters->h.limits) ;
	free(alignment_parameters->h.penalties) ;
	free(alignment_parameters->a.limits) ;
	free(alignment_parameters->a.penalties) ;
	free(alignment_parameters->d.limits) ;
	free(alignment_parameters->d.penalties) ;
	for (int i=0; i<alignment_parameters->num_qualityPlifs; i++)
	{
		free(alignment_parameters->qualityPlifs[i].limits) ;
		free(alignment_parameters->qualityPlifs[i].penalties) ;
	}
	free(alignment_parameters->qualityPlifs) ;
	free(alignment_parameters) ;
	alignment_parameters = NULL ;

	return 0 ;
}


	
int QPalma::read_plif(FILE *fd, struct penalty_struct &plif) {
	char buf[1000] ;

	int narg = fscanf(fd, "%1000s:\t", buf);
	assert(narg==1);
	if (strlen(buf) > 1 && buf[strlen(buf) - 1] == ':')
		buf[strlen(buf) - 1] = 0;
	plif.name = strdup(buf);
	assert(plif.name!=NULL);

	narg = fscanf(fd, "%i\t%i\t\t", &plif.min_len, &plif.max_len);
	//fprintf(stdout, "narg=%i\n", narg) ;
	assert(narg==2);
	//assert(strlen(buf)==0) ;
	plif.transform = T_LINEAR;

	double values[100];
	int values_cnt = 0;
	while (1) {
		narg = fscanf(fd, "%lf,", &values[values_cnt]);
		if (narg != 1)
			break;
		//fprintf(stdout, "%i:%f  ", values_cnt, values[values_cnt]) ;
		values_cnt++;
		assert(values_cnt<100);
	}

	plif.len = values_cnt / 2;
	plif.limits = (double*) malloc(sizeof(double) * plif.len);
	plif.penalties = (double*) malloc(sizeof(double) * plif.len);

	if (plif.limits == NULL || plif.penalties == NULL) {
		fprintf(stderr,	"[read_plif] Could not allocate memory\n");
		free(plif.limits);
		free(plif.penalties);
		return -1;
	}

	for (int i = 0; i < plif.len; i++) {
		plif.limits[i] = values[i];
		plif.penalties[i] = values[plif.len + i];
	}
	plif.cache = NULL;
	plif.use_svm = 0;
	plif.next_pen = NULL;
	return 0;
}

int QPalma::read_matrix(FILE* fd, double *& matrix, char*& name, int dims[2]) {
	char buf[1000];

	int narg = fscanf(fd, "%1000s:\t", buf);
	assert(narg==1);
	if (buf[strlen(buf) - 1] == ':')
		buf[strlen(buf) - 1] = '\0';
	name = strdup(buf);
	assert(name!=NULL);

	narg = fscanf(fd, "%i\t%i\t", &dims[0], &dims[1]);
	assert(narg==2);

	matrix = (double*) malloc(sizeof(double) * dims[0] * dims[1]);

	if (matrix == NULL) {
		fprintf(stderr,
				"[read_matrix] Could not allocate memory (name = %s)\n", name);
		return -1;
	}

	int cnt = 0;
	while (1) {
		narg = fscanf(fd, "%lf,", &matrix[cnt]);
		if (narg != 1)
			break;
		assert(cnt<dims[0]*dims[1]);
		cnt++;
	}
	return 0;
}


int QPalma::init_spliced_align(const char *fname, struct penalty_struct &h,
							   struct penalty_struct &a, struct penalty_struct &d,
							   struct penalty_struct *&qualityPlifs, int &num_qualityPlifs,
							   double*&matchmatrix, int dims[2], int &quality_offset) 
{

	if (_config.NO_QPALMA){
		
		//No scoring function for intron length (h), acceptor sites (a), donor sites (d), edit operations
		init_penalty_struct(h);
		init_penalty_struct(a);
		init_penalty_struct(d);
		qualityPlifs=NULL;
		num_qualityPlifs=0;

		// No quality-based score
		quality_offset=33; //Sanger Fastq format
		
		//Fill matchmatrix from mismatch, gap and match penalties
		//   - A C G T N
		// - x x x x x x
		// A x x x x x x
		// C x x x x x x
		// G x x x x x x
		// T x x x x x x
		// N x x x x x x
		dims[0]=6;
		dims[1]=6;
		
		matchmatrix = (double*) malloc(sizeof(double) * dims[0]*dims[1]);

		if (matchmatrix == NULL) {
			fprintf(stderr,
					"[init_spliced_align] Could not allocate memory\n");
			return -1;
		}
		
		
		for (int i=0;i<dims[0];i++){ //DNA
			for (int j=0;j<dims[1];j++){ //Read

				if (i==0){ // Gap on DNA
					matchmatrix[i*dims[0]+j] = -_config.GAP_SCORE;
					continue;
				}
				if (j == 0){ //Gap on Read
					matchmatrix[i*dims[0]+j] = -_config.GAP_SCORE;
					continue;
				}
				
				if (i==j){ //Match
					matchmatrix[i*dims[0]+j] = -_config.M_SCORE;
					continue;
				}
				
				//Mismatch
				matchmatrix[i*dims[0]+j] = -_config.MM_SCORE;
				
			}
		}
	}
	else
	{
		num_qualityPlifs = 30;
		qualityPlifs = (struct penalty_struct*) malloc(
			sizeof(struct penalty_struct) * num_qualityPlifs);
	
		if (qualityPlifs == NULL) {
			fprintf(stderr, "[init_spliced_align] Could not allocate memory\n");
			return -1;
		}
		char *matrix_name = NULL;
	
		FILE * fd = Util::openFile(fname, "r");
		Util::skip_comment_lines(fd) ;
		int ret = read_plif(fd, h);
		if (ret < 0)
			return ret;

		assert(strcmp(h.name, "h")==0);
		h.max_len = _config.SPLICED_LONGEST_INTRON_LENGTH ;
		for (int i = 0; i < h.len; i++)
			if (h.limits[i] > 2)
				h.limits[i] = 2;

		Util::skip_comment_lines(fd) ;
		ret = read_plif(fd, d);
		if (ret < 0)
			return ret;
		assert(strcmp(d.name, "d")==0);
		d.use_svm = 1;

		Util::skip_comment_lines(fd) ;
		ret = read_plif(fd, a);
		assert(strcmp(a.name, "a")==0);
		if (ret < 0)
			return ret;
		a.use_svm = 1;
		for (int i = 0; i < num_qualityPlifs; i++)
		{
			Util::skip_comment_lines(fd) ;
			read_plif(fd, qualityPlifs[i]);
		}

		Util::skip_comment_lines(fd) ;
		ret = read_matrix(fd, matchmatrix, matrix_name, dims);
		if (ret < 0)
		{
			fprintf(stderr, "init_spliced_align: reading mmatrix failed\n") ;
			return ret;
		}
		assert(strcmp(matrix_name, "mmatrix")==0);
		free(matrix_name);
		matrix_name=NULL ;
	
		double *quality_offset_matrix ;
		int quality_offset_dims[2] ;
		Util::skip_comment_lines(fd) ;
		ret = read_matrix(fd, quality_offset_matrix, matrix_name, quality_offset_dims);
		if (ret < 0)
		{
			fprintf(stderr, "init_spliced_align: reading quality_offset failed\n") ;
			return ret;
		}
		assert(strcmp(matrix_name, "prb_offset")==0);
		free(matrix_name);
		assert(quality_offset_dims[0]==1) ;
		assert(quality_offset_dims[1]==1) ;
		quality_offset=(int)quality_offset_matrix[0] ;
		free(quality_offset_matrix) ;
	}
	
	return 0;
}



int QPalma::get_splicesite_positions(std::string file_template, const char * type, Chromosome const &chr, char strand, int start, int end, float thresh, bool store_pos,
									 std::vector<int> &positions) const
{
	char basename[1000] ;

	// some bounds checking/corrections
	if (start<0)
		start=0 ;
	if (start>(int)chr.length())
		start=chr.length()-1 ;
	
	if (end<0)
		end=0 ;
	if (end>(int)chr.length())
		end=chr.length()-1 ;

	int num=0 ;
	if (_config.NO_SPLICE_PREDICTIONS)
	{
		if (strcmp(type, "acc")==0)
		{
			if (strand=='+')
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.ACC_CONSENSUS.size(); j++)
						if (i>start && chr[i-1] == _config.ACC_CONSENSUS[j][0] && chr[i]==_config.ACC_CONSENSUS[j][1])
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
			else
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.ACC_CONSENSUS_REV.size(); j++)
						if (i>start && chr[i-1] == _config.ACC_CONSENSUS_REV[j][0] && chr[i] == _config.ACC_CONSENSUS_REV[j][1])
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
		}
		if (strcmp(type, "don")==0)
		{
			if (strand=='+')
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.DON_CONSENSUS.size(); j++)
						if (i<end-1 && chr[i]== _config.DON_CONSENSUS[j][0] && chr[i+1]==_config.DON_CONSENSUS[j][1])
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
			else
		    {
				for (int i=start; i<end; i++)
					for (unsigned int j=0; j<_config.DON_CONSENSUS_REV.size(); j++)
						if (i<end-1 && chr[i]== _config.DON_CONSENSUS_REV[j][0] && chr[i+1]==_config.DON_CONSENSUS_REV[j][1])
						{
							num++ ;
							if (store_pos)
								positions.push_back(i) ;
						}
		    } 
		}
	}
	else
	{
		IntervalQuery iq;
		const int num_scores = 1;
		const char *score_names[num_scores] = { "Conf_cum" };
		const int num_intervals = 1 ;
		int interval_matrix[num_intervals * 2];
		interval_matrix[0] = start ;
		interval_matrix[1] = end ;
		
		sprintf(basename, file_template.c_str(), chr.nr() + 1, strand) ;
		
		int ss_size = iq.query(basename, (char**) score_names, num_scores,
							   interval_matrix, num_intervals);
		if (ss_size < 0) 
			return -1;
		
		int* ss_pos = NULL ;
		int* ss_index = NULL ;
		double* ss_score = NULL ;
		try
		{
			//int* ss_pos = (int*) malloc(sizeof(int) * ss_size);
			ss_pos = new int[ss_size] ;
			//int* ss_index = (int*) malloc(sizeof(int) * ss_size);
			ss_index = new int[ss_size] ;
			
			//double* ss_score = (double*) malloc(sizeof(double) * ss_size);
			ss_score = new double[ss_size] ;
		}
		catch (std::bad_alloc&)
		{
			if (ss_pos == NULL || ss_index == NULL || ss_score == NULL) 
			{
				fprintf(stderr, "[map_splice_sites] Could not allocate memory\n");
				delete[] ss_pos;
				delete[] ss_index;
				delete[] ss_score;
				return -1;
			}
		}
		
		
		iq.getResults(ss_pos, ss_index, ss_score);
		iq.cleanup();
		
		for (int i=0; i<ss_size; i++)
		{
			if (ss_score[i]>=thresh)
			{
				num++ ;
				if (store_pos)
				{
					positions.push_back(ss_pos[i]) ;
					//scores.push_back(ss_score[i]) ;
				}
			}
		}
		
		delete[] ss_pos ;
		delete[] ss_index ;
		delete[] ss_score ;
	}
	
	return num ;
}

int QPalma::get_string_from_region(Chromosome const &chrN, region_t *region, std::string &str) const
{

	assert(region->end > region->start);

	int32_t length = region->end - region->start;

	str.assign("") ;
	for (int i=0; i<length; i++)
		str.push_back(chrN[region->start+i]) ;
	//str.assign(CHR_SEQ[chrN] + region->start, length);

	// add (up to) three N letters on each side
	if (str.size()>0)
	{
		str[0]='N' ;
		str[str.size()-1]='N' ;
	}
	if (str.size()>1)
	{
		str[1]='N' ;
		str[str.size()-2]='N' ;
	}
	if (str.size()>2)
	{
		str[2]='N' ;
		str[str.size()-3]='N' ;
	}
	return 0 ;
}

void QPalma::Result::add_buffer_to_region(int ori, Chromosome const &chrN, int32_t nregion) {
	region_t *region = regions[ori][chrN.nr()][nregion];

	if (region->start <= _config.SPLICED_CLUSTER_TOLERANCE)
		region->start = 1;
	else
		region->start -= _config.SPLICED_CLUSTER_TOLERANCE;

	if (region->end + _config.SPLICED_CLUSTER_TOLERANCE >= (int)chrN.length())
		region->end = chrN.length();// - 1; chr length because end bound is excluded
	else
		region->end += _config.SPLICED_CLUSTER_TOLERANCE;
}

/** performs a quicksort on an array output of length size
 * it is sorted from in ascending (for type T) */
void QPalma::qsort(region_t** output, int size) const {
	if (size == 2) {
		if (output[0]->start > output[1]->start) {
			region_t *c = output[0];
			output[0] = output[1];
			output[1] = c;
		}
		return;
	}
	//T split=output[random(0,size-1)];
	region_t *split = output[size / 2];

	int32_t left = 0;
	int32_t right = size - 1;

	while (left <= right) {
		while (output[left]->start < split->start)
			left++;
		while (output[right]->start > split->start)
			right--;

		if (left <= right) {
			region_t *c = output[left];
			output[left] = output[right];
			output[right] = c;
			left++;
			right--;
		}
	}

	if (right + 1 > 1)
		qsort(output, right + 1);

	if (size - left > 1)
		qsort(&output[left], size - left);
}


void QPalma::print_map(Read const &read, bool* read_map, const char *name)
{
	fprintf(stdout, "# read_map (%s): ", name);
	for (size_t i = 0; i < read.length(); i++)
		if (read_map[i])
			fprintf(stdout, "1");
		else
			fprintf(stdout, "0");
	fprintf(stdout, "\n");
}


////////////////////////////////////////////////////////////////////////////////
// filter


int QPalma::get_num_splicesites(std::string file_template, const char* type, Chromosome const &chr, char strand, int start, int end, float thresh) const
{
	std::vector<int> positions ;
	return get_splicesite_positions(file_template, type, chr, strand, start, end, thresh, false, positions) ;
}


int QPalma::qpalma_filter(Result &result, struct alignment_t *ali, int num_N,JunctionMap &annotatedjunctions) const
{
	assert(ali->exons.size()<=2) ;
	static bool use_ss_for_filtering = true ;
	
	Chromosome const &chr = *ali -> chromosome ;
	int start = ali -> exons[0] ;
	int end = ali -> exons[1] ;
	
	unsigned int num_gaps = ali -> num_gaps ;
	unsigned int num_matches = ali -> num_matches ;

	if (num_gaps > _config.FILTER_BY_MAX_GAPS || result._read.length()-num_matches > _config.FILTER_BY_MAX_MISMATCHES+num_N)
	{
		if (verbosity>=1)
			fprintf(stdout, "filter decides YES: num_gaps=%i, num_mismatches=%i, num_N=%i\n", num_gaps, result._read.length()-num_matches, num_N) ;
		
		result.qpalma_filter_reason=2 ; // cheap positive filter
			
		return 1 ;
	}

	float thresh_acc = _config.FILTER_BY_SPLICE_SITES_THRESH_ACC ;
	float thresh_don = _config.FILTER_BY_SPLICE_SITES_THRESH_DON ;
	int region = _config.FILTER_BY_SPLICE_SITES_REGION ;

	if (region==-1 || !_config.FILTER_BY_SPLICE_SITES)
	{
		result.qpalma_filter_reason=0 ; // cheap negative filter
		return 0 ;
	}

	if ( result._read.length()-num_matches < _config.FILTER_BY_SPLICE_SITES_EDIT_MIN+num_N )
	{
		result.qpalma_filter_reason=0 ; // cheap negative filter
		return 0 ;
	}
	
	if (use_ss_for_filtering)
		try
		{
			int num_ss = 0;
			
			if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0){
				

				num_ss+= get_num_splicesites(_config.ACC_FILES, "acc", chr, '+', start-region, end+region, thresh_acc) + 
				get_num_splicesites(_config.ACC_FILES, "acc", chr, '-', start-region, end+region, thresh_acc) +
				get_num_splicesites(_config.DON_FILES, "don", chr, '+', start-region, end+region, thresh_don) +
				get_num_splicesites(_config.DON_FILES, "don", chr, '-', start-region, end+region, thresh_don) ;
			}
			if (_config.SCORE_ANNOTATED_SPLICE_SITES){
				std::vector<int> positions;
				int chrN = chr.nr();
				
				get_annotated_splice_positions(positions,annotatedjunctions,"acc",start-region, end+region,chrN,'+');
				get_annotated_splice_positions(positions,annotatedjunctions,"don",start-region, end+region,chrN,'+');
				get_annotated_splice_positions(positions,annotatedjunctions,"acc",start-region, end+region,chrN,'-');
				get_annotated_splice_positions(positions,annotatedjunctions,"don",start-region, end+region,chrN,'-');
				num_ss += positions.size();
				//fprintf(stdout,"Positions %i\n",positions.size());
				
				positions.clear();
				

			}
			
			
			if ( num_ss>0 )
			{
				if (verbosity>=1)
					fprintf(stdout, "filter decides YES: num_ss=%i, gaps=%i, num_mismatches=%i, num_N=%i\n", num_ss, num_gaps, result._read.length()-num_matches, num_N) ;
				
				result.qpalma_filter_reason=3 ; // positive splice site filter
				
				return 1 ;
			}
		}
		catch (IntervalQueryException & e)
		{
			fprintf(stdout, "Warning: do not use splice site files for triggering qpalma alignments\n") ;
			use_ss_for_filtering=false ;
		}
	
	if (verbosity>=1)
		fprintf(stdout, "filter decides NO: num_gaps=%i, num_mismatches=%i, num_N=%i\n", num_gaps, result._read.length()-num_matches, num_N) ;
	
	result.qpalma_filter_reason=1 ; // negative splice site filter

	return 0 ;
}




////////////////////////////////////////////////////////////////////////////////
// alignment

/** find long regions included in the set of current regions */
void QPalma::recover_long_regions(Read const &read, std::vector<region_t*> &long_regions_output, std::vector<region_t*> long_regions, std::vector<region_t*> current_regions) const
{

	// Sort long_regions by start position
	region_t ** arr = NULL ;
	size_t nbr_long_regions=long_regions.size();
	try 
    {
		arr = new region_t*[nbr_long_regions] ;			
    }
	catch (std::bad_alloc&)
    {
		fprintf(stderr, "[capture_hits] ERROR Could not allocate memory (_read.id() = %s)\n",
				read.id());
    }
	for (int i = 0; i < (int)nbr_long_regions; i++)
		arr[i] = long_regions[i];  
	qsort(arr, nbr_long_regions);
 
	unsigned int reg=0;
	unsigned int lr=0;
  
	//fprintf(stdout,"**start** %i-%i\n",current_regions.size(),nbr_long_regions);
	while (reg < current_regions.size() && lr < nbr_long_regions){
		//fprintf(stdout,"current region: %i-%i and current long region: %i-%i\n", current_regions[reg]->start,current_regions[reg]->end,arr[lr]->start,arr[lr]->end);
    
		if (arr[lr]->start >= current_regions[reg]->start && arr[lr]->end <= current_regions[reg]->end){
			//print_map(arr[lr]->read_map,"array");
			long_regions_output.push_back(arr[lr]);
			//print_map(long_regions_output[long_regions_output.size()-1]->read_map,"copy ");
			//fprintf(stdout,"Add one long region as starting point\n");
			lr++;
		}
		else if (arr[lr]->start < current_regions[reg]->start)
			lr++;
		else // current long region coordinates after current region
			reg++;
	}
	delete[] arr;
	//fprintf(stdout,"**stop** %i-%i\n",reg,lr);
}

/** Gives the relative position on dna sequence to align  */
int QPalma::convert_dna_position(int real_position, size_t* cum_length, const std::vector<region_t *> &current_regions) const
{
	for (size_t j = 0; j < current_regions.size(); j++)
		if (real_position >= current_regions[j]->start && real_position< current_regions[j]->end)
			return real_position - current_regions[j]->start + cum_length[j];

	//   for (int i=0;i<positions.size();i++){
//       if (positions[i]==real_position)
// 	return i;
//     }
  
	return -1;
}

int QPalma::get_first_read_map(Read const &read, bool* read_map) const
{

	for(unsigned int i=0; i < read.length(); i++)
		if (read_map[i])
			return i;
  
	return -1;
}

// GenomeMapper's handling of spliced reads is very simplistic. It just finds
// exact 12-mers and returns these as spliced hits. It does _not_ do an
// alignment of any kind. The following routine takes a pair of such
// spliced hits and aligns them against the region of the genome where the
// hits were reported.


void QPalma::print_hit(HIT *hit) {
	fprintf(stdout, "# hit= {start=%i, end=%i, readpos=%i}\n", hit->start,
			hit->end, hit->readpos);
}

void QPalma::print_region(region_t *region, const char * bla) 
{
	fprintf(stdout, "# region %s = {start=%i, end=%i, ori=%c, erased=%i}\n", bla,
			region->start, region->end, region->orientation, region->erased);
}


void QPalma::Result::delete_regions()
{
	for (int ori = 0; ori < 2; ori++)
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++)
			for (int32_t nregion = 0; nregion < (int)regions[ori][chrN].size(); nregion++) 
			{
				delete[] regions[ori][chrN][nregion]->read_map;
				regions[ori][chrN][nregion]->read_map = NULL;
				delete regions[ori][chrN][nregion];
			}
}

void QPalma::delete_long_regions(std::vector<std::vector<region_t *> > *long_regions) const
{ 
	//long_regions need to be deleted: deep copy of initial regions into long_regions
	for (int ori = 0; ori < 2; ori++)
		for (int32_t chrN = 0; chrN < (int)long_regions[ori].size(); chrN++)
			for (int32_t nregion = 0; nregion < (int)long_regions[ori][chrN].size(); nregion++) 
			{
				delete[] long_regions[ori][chrN][nregion]->read_map;
				long_regions[ori][chrN][nregion]->read_map = NULL;
				delete long_regions[ori][chrN][nregion];
			}
}


int QPalma::capture_hits(Hits &hits, Result &result, bool const non_consensus_search,JunctionMap &annotatedjunctions) const
{
	Read const &read(hits.getRead());

	clock_t start_time = clock();
  

	// clean up data generated for the previous read

	result.cleanup();

	HIT const *hit;
	int32_t num_hits = 0; // TODO debugging only
	int32_t num_hits_dropped = 0; // TODO debugging only
  
	// Examine all hits and construct a list of region where these hits map.
	// regions is a list of region clusters per chromosome, sorted by start position
  
	result.long_regions[0] = result.regions[0];
	result.long_regions[1] = result.regions[1];
	std::vector<std::vector<region_t *> > * const regions = result.regions;
	std::vector<std::vector<region_t *> > * const long_regions = result.long_regions;
  

	//TODO: Real length of a hit is i-1
	for (int32_t i = read.length(); i >= _config.SPLICED_HIT_MIN_LENGTH_SHORT; i--) {

		
		hit = hits.HIT_LISTS_OPERATOR[i];
    
		while (hit != NULL) 
		{
			//bool captured = false;
			
			num_hits++;
			//TODO: Real length of a hit is i-1
			bool consider_as_long_hit = (i >= _config.SPLICED_HIT_MIN_LENGTH_LONG);
			if (consider_as_long_hit && num_hits >= MAX_NUM_LONG_HITS) {
				consider_as_long_hit = false;
				if (verbosity >= 2)
					fprintf(stdout, "# ignoring long hits from now on\n");
			}
	

			if (!consider_as_long_hit) {
				// first check whether it is close enough to a long enough hit
				bool found = false;
				for (int32_t nregion = 0; nregion < (int)long_regions[ori_map(hit->orientation)][hit->chromosome->nr()].size(); nregion++) {
					//if (!long_regions[hit->chromosome][nregion]->islong)
					//	continue ;
					int32_t rs = long_regions[ori_map(hit->orientation)][hit->chromosome->nr()][nregion]->start - _config.SPLICED_LONGEST_INTRON_LENGTH;
					if (rs < 0)
						rs = 0;
					int32_t re = long_regions[ori_map(hit->orientation)][hit->chromosome->nr()][nregion]->end + _config.SPLICED_LONGEST_INTRON_LENGTH;
					if (((int)hit->start >= rs) && ((int)hit->end <= re)) 
					{
						found = true;
						break;
					}
				}
				if (!found) {
					num_hits_dropped++;
					//	    fprintf(stdout, "dropped %i-%i (%i)\n", hit->start, hit->end,hit->readpos) ;
					hit = hit->next;
					continue;
				}
			}
			if (verbosity >= 2){
				fprintf(stdout, "# captured %i-%i  (%i,%i) %c %s real length: %i mism-gaps(%i-%i)\n", hit->start,
						hit->end, hit->readpos, hit->end - hit->start, 
						hit->orientation,consider_as_long_hit?"long":"notlong",i,hit->mismatches,hit->gaps);
// 	  if (consider_as_long_hit){
// 	    for(int n=hit->start;n<hit->end;n++)
// 	      fprintf(stdout,"%c",CHR_SEQ(hit->chromosome,n));
// 	    fprintf(stdout,"\n");
// 	    for(int n=0;n<hit->end-hit->start;n++){
// 	      fprintf(stdout,"%c",READ[hit->readpos+n]);
// 	    }
// 	    fprintf(stdout,"\n");
// 	    fprintf(stdout,"%s\n",READ);
// 	  }
			}
			//if (regions[ori_map(hit->orientation)][hit->chromosome].empty() && consider_as_long_hit)
			{
				// Create first region for this hit
				region_t *new_region = NULL;
				try {
					new_region = new region_t();
					new_region->read_map = new bool[read.length()];
				} catch (std::bad_alloc&) 
				{
					fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
					result.delete_regions();
					delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
					return -1;
				}  
	  
				new_region->start = hit->start;
				new_region->end = hit->end;
				new_region->from_map = false ;
				for (size_t ii = 0; ii < read.length(); ii++)
					new_region->read_map[ii] = false;
				assert(hit->end >= hit->start) ;

				
				for (size_t ii = 0; ii < hit->end - hit->start && hit->readpos + ii < read.length(); ii++)
					new_region->read_map[hit->readpos + ii] = true;
				//print_map(new_region->read_map, "init") ;
				//new_region->strand = hit->orientation ;
				//new_region->chromosome = hit->chromosome ;
	  
				if (consider_as_long_hit){
					region_t *new_lregion=NULL;
					try {
						new_lregion = new region_t();
						new_lregion->read_map = new bool[read.length()];
					} catch (std::bad_alloc&) 
					{
						fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
						result.delete_regions();
						delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
						return -1;
					}  
	  
					new_lregion->start = hit->start;
					new_lregion->end = hit->end;
					new_lregion->from_map = false ;
					for (size_t ii = 0; ii < read.length(); ii++)
						new_lregion->read_map[ii] = false;
					for (size_t ii = 0; ii < hit->end - hit->start && hit->readpos + ii< read.length(); ii++)
						new_lregion->read_map[hit->readpos + ii] = true;
					long_regions[ori_map(hit->orientation)][hit->chromosome->nr()].push_back(new_lregion);
				}
				regions[ori_map(hit->orientation)][hit->chromosome->nr()].push_back(new_region);
	  

				// 
				hit = hit->next;
				continue;
			}
	
			/*for (int32_t tmp = 0; tmp < regions[hit->chromosome].size(); tmp++)
			  fprintf(stdout, "\t%d - %d  [%i]\n", regions[hit->chromosome][tmp]->start,
			  regions[hit->chromosome][tmp]->end, regions[hit->chromosome][tmp]->islong);
			  fprintf(stdout, "\n") ;*/
	
			hit = hit->next;
		}
	}
	if (verbosity >= 1)
		fprintf(stdout,	"# [capture_hits] Captured %d hits, dropped %i hits for read %s\n",
				num_hits - num_hits_dropped, num_hits_dropped, read.id());
  
	//find short hits in close vicinity
	//find_short_hits();
  
	if (_config.QPALMA_USE_MAP || _config.QPALMA_USE_SPLICE_SITES)
    {
		// add regions from reporting map
		//const int take_report_map = MASK_MAPPED_READ_BEST | MASK_MAPPED_READ | MASK_SPLICED_READ_BEST | MASK_SPLICED_READ | MASK_MAPPED_REGION ;
		//const int take_report_map = 0 ;
      
		//const int take_report_map = MASK_MAPPED_READ_BEST | MASK_SPLICED_READ_BEST ;
		const int take_report_map = MASK_MAPPED_READ_BEST | MASK_SPLICED_READ_BEST ; 
      
		int added_map_regions = 0 ; 
		int added_map_total_len = 0 ;
      
		for (int ori = 0; ori < 2; ori++)
		{
			for (int32_t chrN = 0; chrN < (int)long_regions[ori].size(); chrN++)
			{
				Chromosome const &chromosome = genome->chromosome(chrN);
				if (long_regions[ori][chrN].size() == 0)
					continue;
				for (int i=0; i<(int)long_regions[ori][chrN].size(); i++)
				{
					int start = long_regions[ori][chrN][i]->start - _config.SPLICED_LONGEST_INTRON_LENGTH ;
					if (start<0)
						start=0 ;
					int end = long_regions[ori][chrN][i]->end + _config.SPLICED_LONGEST_INTRON_LENGTH ;
					if (end>(int)chromosome.length())
						end=chromosome.length();
					int midpoint = (end+start)/2 ;
		  
					
					if (_config.QPALMA_USE_SPLICE_SITES && !_config.NO_SPLICE_PREDICTIONS)

					{
						std::vector<int> positions;
						int num_acc=0, num_don=0 ;
						
						// find all acceptor splice sites downstream and donor splice site upstream
						if (ori==0)
						{
							if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0){
								num_acc = get_splicesite_positions(_config.ACC_FILES, "acc", chromosome, '+', midpoint, end, _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC, true, positions) ;
								num_don = get_splicesite_positions(_config.DON_FILES, "don", chromosome, '+', start, midpoint, _config.QPALMA_USE_SPLICE_SITES_THRESH_DON, true, positions) ;
							}
							
							if (_config.SCORE_ANNOTATED_SPLICE_SITES){
								//fprintf(stdout,"Annotated splice site for pseudo chromosome seq\n");
								
								get_annotated_splice_positions(positions,annotatedjunctions,"acc",midpoint,end,chrN,'+');
								get_annotated_splice_positions(positions,annotatedjunctions,"don",start,midpoint,chrN,'+');
							}
						}
						
						else
						{
							if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0){
								num_don = get_splicesite_positions(_config.DON_FILES, "don", chromosome, '-', midpoint, end, _config.QPALMA_USE_SPLICE_SITES_THRESH_DON, true, positions) ;
								num_acc = get_splicesite_positions(_config.ACC_FILES, "acc", chromosome, '-', start, midpoint, _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC, true, positions) ;
								}
							if (_config.SCORE_ANNOTATED_SPLICE_SITES){
								get_annotated_splice_positions(positions,annotatedjunctions,"acc",midpoint,end,chrN,'-');
								get_annotated_splice_positions(positions,annotatedjunctions,"don",start,midpoint,chrN,'-');
							}	
							
							
						}
						
						if (verbosity>=1)
							fprintf(stdout, "splice site map: chr=%i, ori=%i: num_acc=%i\t num_don=%i\n", chrN, ori, num_acc, num_don) ;

						for (size_t jj=0; jj<positions.size(); jj++)
						{
							region_t *new_region = NULL;
							try {
								new_region = new region_t();
								new_region->read_map = new bool[read.length()];
							} catch (std::bad_alloc&) {
								fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
								result.delete_regions();
								delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
								return -1;
							}
			  
							new_region->start = positions[jj] - _config.INDEX_DEPTH ;
							new_region->end = positions[jj] + _config.INDEX_DEPTH ;
							for (size_t ii = 0; ii < read.length(); ii++)
								new_region->read_map[ii] = false;

							//new_region->strand = (ori == 0) ? '+' : '-' ;
							//new_region->chromosome = chrN ;
							new_region->from_map = true ;
							
							regions[ori][chrN].push_back(new_region);
							assert(new_region->end > new_region->start) ;
			  
							//fprintf(stdout, "ss region %i-%i (seed %i-%i)\n", new_region->start, new_region->end, start, end) ;
			  
							added_map_regions++ ;
							added_map_total_len += new_region->end - new_region->start; 
						}
					}
		  
					if (_config.QPALMA_USE_MAP)
					{
						size_t region_ptr = regions[ori][chrN].size() ;
		      
						int region_start = -1 ;
						int region_end = -1 ;
						for (int p=start; p<end; p+=Config::QPALMA_USE_MAP_WINDOW)//p+=1)
						{
							int val = genomemaps->CHR_MAP(chromosome,p) ;
			  
							if ((val & take_report_map)!=0 && end-p>Config::QPALMA_USE_MAP_WINDOW)//p<end-1)
							{
								if (region_start==-1){
									region_start=p;
								}
								region_end=p;
							}
							else
							{
								if (region_start!=-1)
								{ // region finished
				  
									region_t *new_region = NULL;
									try {
										new_region = new region_t();
										new_region->read_map = new bool[read.length()];
									} catch (std::bad_alloc&) {
										fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n");
										result.delete_regions();
										delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
										return -1;
									}
				  
									new_region->start = region_start - _config.INDEX_DEPTH ;
									if (end-p<=Config::QPALMA_USE_MAP_WINDOW)
										new_region->end = end-1 + _config.INDEX_DEPTH ;
									else
										new_region->end = region_end + _config.INDEX_DEPTH ;
									for (size_t ii = 0; ii < read.length(); ii++)
										new_region->read_map[ii] = false;
									//new_region->strand = (ori == 0) ? '+' : '-' ;
									//new_region->chromosome = chrN ;
									new_region->from_map = true ;
				  
									regions[ori][chrN].push_back(new_region);
									assert(new_region->end > new_region->start) ;
								}
								region_start = -1 ;
								region_end = -1 ;
							}
						}
						
						
						// look at the regions downstream. only keep those which lead to an alignment string shorter than _config.QPALMA_USE_MAP_MAX_SIZE
						int right_covered_region=0 ;
						int right_removed = 0 ;
						for (size_t nregion=region_ptr; nregion<regions[ori][chrN].size(); nregion++)
						{
							assert(regions[ori][chrN][nregion]->from_map) ;
							if (regions[ori][chrN][nregion]->start > midpoint)
							{
								//fprintf(stdout, "region %i: %i - %i\n", (int)nregion, regions[ori][chrN][nregion]->start, regions[ori][chrN][nregion]->end) ;
								if (right_covered_region>_config.QPALMA_USE_MAP_MAX_SIZE)
								{
									regions[ori][chrN][nregion]->erased = true ;
									delete[] regions[ori][chrN][nregion]->read_map ;
									regions[ori][chrN][nregion]->read_map=NULL ;
									right_removed++ ;
									//fprintf(stdout, "dropped\n") ;
								}
								else
								{
									added_map_regions++ ;
									added_map_total_len += regions[ori][chrN][nregion]->end - regions[ori][chrN][nregion]->start;
									//fprintf(stdout, "accepted\n") ;
								}
								right_covered_region+= regions[ori][chrN][nregion]->end - regions[ori][chrN][nregion]->start;
							}
						}
						// look at the regions upstream. only keep those which lead to an alignment string shorter than _config.QPALMA_USE_MAP_MAX_SIZE
						int left_covered_region=0 ;
						int left_removed=0 ;
						for (size_t nregion=regions[ori][chrN].size()-1; nregion>=region_ptr; nregion--)
						{
							assert(regions[ori][chrN][nregion]->from_map) ;
							if (regions[ori][chrN][nregion]->start < midpoint)
							{
								//fprintf(stdout, "region %i: %i - %i\n", (int)nregion, regions[ori][chrN][nregion]->start, regions[ori][chrN][nregion]->end) ;
								if (left_covered_region > _config.QPALMA_USE_MAP_MAX_SIZE)
								{
									regions[ori][chrN][nregion]->erased = true ;
									delete[] regions[ori][chrN][nregion]->read_map ;
									regions[ori][chrN][nregion]->read_map=NULL ;
									left_removed++ ;
									//fprintf(stdout, "dropped\n") ;
								}
								else
								{
									added_map_regions++ ;
									added_map_total_len += regions[ori][chrN][nregion]->end - regions[ori][chrN][nregion]->start;
									//fprintf(stdout, "accepted\n") ;
								}
								left_covered_region += (regions[ori][chrN][nregion]->end - regions[ori][chrN][nregion]->start);
							}
						}
						if (verbosity>=2)
							fprintf(stdout, "cover: %i, %i\tregions removed: %i, %i\n", left_covered_region, right_covered_region, left_removed, right_removed) ;
					}
				}
			}
		}
		if (verbosity>=1)
			fprintf(stdout, "# added %i regions from reporting map (total len=%i, average len=%i)\n", added_map_regions, added_map_total_len, added_map_total_len/(added_map_regions+1)) ;
    }
  
	// sort regions by starting position
	for (int ori = 0; ori < 2; ori++)
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++) {
			if (regions[ori][chrN].size() == 0)
				continue;
      
			region_t ** arr = NULL ;
			size_t nbr_regions=regions[ori][chrN].size();
			try 
			{
				arr = new region_t*[nbr_regions] ;			
			}
			catch (std::bad_alloc&)
			{
				fprintf(stderr, "[capture_hits] ERROR Could not allocate memory (_read.id() = %s)\n",
						read.id());
				return -1;
			}
      
			for (int i = 0; i < (int)nbr_regions; i++)
				arr[i] = regions[ori][chrN][i];
      
			qsort(arr, nbr_regions);
      
			for (int i = 0; i < (int)nbr_regions; i++)
				regions[ori][chrN][i] = arr[i];
      
			for (int i = 0; i + 1 < (int)nbr_regions; i++)
				assert(regions[ori][chrN][i]->start<=regions[ori][chrN][i+1]->start);
      
			//for (int i=0; i<nbr_regions; i++)
			//	arr[i]=NULL ;
      
			delete[] arr;
		}

	// Surround the regions with a buffer of size _config.SPLICED_CLUSTER_TOLERANCE.
	// If any of them overlap before/after extension, merge them to
	// avoid duplicating subsequences.
	//
    
    for (int ori = 0; ori < 2; ori++){
		for (int32_t chrN = 0; chrN < (int)regions[ori].size(); chrN++)
		{
			Chromosome const &chromosome = genome->chromosome(chrN);
			size_t nbr_regions=regions[ori][chrN].size();
			if (nbr_regions == 0)
				continue;
			if (nbr_regions == 1) 
			{
				// Nothing to merge, just extend the one existing region to include a buffer.
				result.add_buffer_to_region(ori, chromosome, 0);
				continue;
			}
	
			for (int i = 0; i < (int)nbr_regions; i++) 
				result.add_buffer_to_region(ori, chromosome, i);
	
			for (int nregion = 0; nregion < (int)nbr_regions - 1; nregion++) 
			{
				if (regions[ori][chrN][nregion]->erased)
					continue;
				size_t next = nregion + 1;
				while (next < nbr_regions && regions[ori][chrN][next]->erased)
					next++;
				if (next >= nbr_regions || regions[ori][chrN][next]->erased)
					continue;
	    
				if ((int32_t) ((regions[ori][chrN][nregion])->end) >= (int32_t) ((regions[ori][chrN][next])->start)){
					// regions[ori] are overlapping or adjacent. Merge them into one.
					if ((regions[ori][chrN][nregion])->end < (regions[ori][chrN][next])->end)
						(regions[ori][chrN][nregion])->end = (regions[ori][chrN][next])->end;
		
					// merge the two read_maps
					for (size_t i = 0; i < read.length(); i++)
						regions[ori][chrN][nregion]->read_map[i] = regions[ori][chrN][nregion]->read_map[i]
							|| regions[ori][chrN][next]->read_map[i];
		
					// Get rid of the extra region.
					regions[ori][chrN][next]->erased = true;
					delete[] regions[ori][chrN][next]->read_map; 
					regions[ori][chrN][next]->read_map = NULL;
		
					// consider this as a map-region, only if both were generated from a map
					regions[ori][chrN][nregion]->from_map = regions[ori][chrN][nregion]->from_map && regions[ori][chrN][next]->from_map ;
		
					// make sure this item is looked at again to merge it with the next item if necessary
					nregion-- ;
				}
			}
		}

		/*
		  for (int nregion=0; nregion<regions[ori][chrN].size(); nregion++)
		  {
		  if (regions[ori][chrN][nregion]->erased)
		  continue ;
		  print_map(regions[ori][chrN][nregion]->read_map, "region") ;
		  }
		*/

	}
  
	//TODO: mt region1_time += clock() - start_time;
  
	
	/*printf("Clusters after consolidation: \n");
	  for (int ori=0; ori<1; ori++)
	  for (int32_t chrN = 0; chrN < regions[ori].size(); chrN++)
	  for (int32_t tmp = 0; tmp < regions[ori][chrN].size(); tmp++)
	  {
	  if (regions[ori][chrN][tmp]->erased)
	  continue ;
	  printf("\t%d - %d\n", regions[ori][chrN][tmp]->start, regions[ori][chrN][tmp]->end);
	  assert(regions[ori][chrN][tmp]->start<regions[ori][chrN][tmp]->end) ;
	  if (tmp>0)
	  assert(regions[ori][chrN][tmp-1]->end<regions[ori][chrN][tmp]->start) ;
	  }
	  fprintf(stdout, "_config.SPLICED_CLUSTER_TOLERANCE=%i\n", _config.SPLICED_CLUSTER_TOLERANCE) ;*/
  
	int num_merged = 0;
	for (int ori = 0; ori < 2; ori++)
		for (size_t chrN = 0; chrN < regions[ori].size(); chrN++)
			for (size_t tmp = 0; tmp < regions[ori][chrN].size(); tmp++)
			{
				assert(regions[ori][chrN][tmp]->end>regions[ori][chrN][tmp]->start) ;
				if (!regions[ori][chrN][tmp]->erased)
					num_merged++;
			}
	if (verbosity >= 2)
		fprintf(stdout, "# Merged hit list has %i items\n", num_merged);
	_stats.qpalma_region_align_time += clock() - start_time;
	return 0;
}

int QPalma::capture_hits_2(Hits &hits, Result &result, bool non_consensus_search, JunctionMap &annotatedjunctions) const 
{
	// Dominik marker
	// std::vector<std::vector<region_t *> > long_regions_other[2] ;
	// if (paired_read_processing)
	//      	long_regions_other = get_long_regions_from_other_end(long_regions, semaphore_whatever_is_necessary) ;

	Read const &read(hits.getRead());
	clock_t start_time = clock();

	std::vector<std::vector<region_t *> > *regions = result.regions;
	std::vector<std::vector<region_t *> > *long_regions = result.long_regions;

	// Construct a pseudo genomic sequence from the regions by stringing
	// them together. A new sequence is started when two adjacent regions
	// are farther apart than the maximum assumed intron length.
	// Also, we build up a map of positions on the generated strings to
	// their corresponding positions in the genomic sequence. This is
	// so we can map the results of the alignment back onto the actual
	// genomic sequence later.

	int num_alignments_reported = 0 ;
	

  
	std::string read_seq[2];
	read_seq[0] = std::string(read.data(), read.length());
	read_seq[1] = reverse(complement(read_seq[0]));
	if (verbosity >= 3) {
		fprintf(stdout, "# read[0]: %s\n", read_seq[0].c_str());
		fprintf(stdout, "# read[1]: %s\n", read_seq[1].c_str());
	}
  
	std::string read_quality[2];
	read_quality[0] = std::string(read.quality(0), read.length());
	read_quality[1] = reverse(read_quality[0]);
	if (verbosity >= 3)
		fprintf(stdout, "# readqual[0]: %s\n", read_quality[0].c_str());
  
	bool *read_map = NULL;
	try 
    {
		read_map = new bool[read.length()];
    } 
	catch (std::bad_alloc&)
    {
		fprintf(stderr, "[capture_hits] allocating memory for read_map failed\n") ;
		result.delete_regions() ;
		delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
		return -1 ;
    }

	for (int ori = 0; ori < 2; ori++){

		hits.topAlignments().init_top_alignment_indice();
	
		for (size_t chrN = 0; chrN < regions[ori].size(); chrN++) 
		{
			Chromosome const &chr = genome->chromosome(chrN);
			std::vector<region_t *> current_regions;
			current_regions.clear();
			std::vector<int> current_positions;
			current_positions.clear();
			if (regions[ori][chrN].size() == 0)
				continue;
			int start_region = -1;
			for (size_t nregion = 0; nregion < regions[ori][chrN].size(); nregion++)
				if (!regions[ori][chrN][nregion]->erased) 
				{
					start_region = nregion;
					break;
				}
			if (start_region == -1)
				continue;
		  
			std::string str;
			{
				int ret = get_string_from_region(chr, regions[ori][chrN][start_region], str);
				if (ret < 0)
				{
					result.delete_regions();
					delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
					return ret;
				}
			}
		  
			current_regions.push_back(regions[ori][chrN][start_region]);
			std::string current_seq = str;
		  
			assert(regions[ori][chrN][start_region]->end >= regions[ori][chrN][start_region]->start) ;
			for (int p = 0; p < regions[ori][chrN][start_region]->end - regions[ori][chrN][start_region]->start; p++) 
			{
				current_positions.push_back(regions[ori][chrN][start_region]->start + p);
			}
		  
			// if (paired_read_processing) 
			//      	"do some clever processing of the region list based on long_regions_other" ;


			// initialize read_map
			for (size_t i = 0; i < read.length(); i++)
				read_map[i] = regions[ori][chrN][start_region]->read_map[i];
		  
			for (size_t nregion = start_region + 1; nregion < regions[ori][chrN].size(); nregion++) 
			{
				if (regions[ori][chrN][nregion]->erased)
					continue;
				{
					int ret = get_string_from_region(chr, regions[ori][chrN][nregion], str);
					if (ret < 0)
					{
						result.delete_regions();
						delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
						return ret;
					}
				}
			  
				if (regions[ori][chrN][nregion]->start - regions[ori][chrN][nregion - 1]->end > _config.SPLICED_LONGEST_INTRON_LENGTH) 
				{
					// Regions are too far apart to contain parts of one spliced
					// hit. Start a new sequence which will be treated in a
					// separate alignment.
				  
					int num_read_map = 0;
					for (size_t i = 0; i < read.length(); i++)
						if (read_map[i])
							num_read_map++;
				  
					if (verbosity >= 2)
						fprintf(stdout, "# region list covering %i bases\n",
								num_read_map);
					if (num_read_map >= _config.SPLICED_HIT_MIN_LENGTH_COMB)// && current_regions.size()>=2)
					{
					  
						//Recover long regions (starting points in sequences) for current regions to align
						std::vector<region_t*> corres_long_regions;
						corres_long_regions.clear();
						recover_long_regions(read, corres_long_regions, long_regions[ori][chrN], current_regions);
						assert(corres_long_regions.size()>0); // at least one long region as support of alignment
						assert(corres_long_regions[0]->read_map!=NULL);
					  
						//Take the first long region  to start alignment
						int hit_read_position = get_first_read_map(read, corres_long_regions[0]->read_map);
						int hit_len= corres_long_regions[0]->end-corres_long_regions[0]->start;
						if(ori==1){
							hit_read_position = read.length()-hit_len-hit_read_position+1;
						}
						assert (hit_read_position>=0 && hit_len >0);
						// fprintf(stdout,"read id %s curr len %i\n",read.id(), current_positions.size());
						// fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",hit_read_position, corres_long_regions[0]->start, hit_len);					  
						// fprintf(stdout,	"# Number of current regions %i\n",(int)current_regions.size());
						
						bool isunspliced ;

						int transcription_direction = get_transcription_direction(_config.STRAND,ori);
						
						if (transcription_direction >=0)
						{
							int ret = perform_alignment_starter(result, hits, read_seq[ori], read_quality[ori], current_seq, current_regions,
																current_positions, chr, '+', ori, hit_read_position,
																corres_long_regions[0]->start, hit_len, non_consensus_search, num_alignments_reported,false,annotatedjunctions);
							if (ret < 0)
							{
								result.delete_regions();
								delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
								return ret;
							}
							isunspliced = ret;
						}
					  
						//if (!isunspliced) 
						if (transcription_direction <=0)
						{
							//fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",_read.lenght()-(hit_read_position+hit_len),
							//      corres_long_regions[0]->end, hit_len);					  
							//fprintf(stdout,	"# Number of current regions %i\n",(int)current_regions.size());					  
							int ret = perform_alignment_starter(result, hits, read_seq[1 - ori], read_quality[1 - ori], current_seq, current_regions, 
																current_positions, chr, '-', 1-ori, read.length()-(hit_read_position+hit_len),
																corres_long_regions[0]->end-1, hit_len, non_consensus_search, num_alignments_reported,false,annotatedjunctions);//end nucleotide in dna not included
							if (ret < 0)
							{
								result.delete_regions();
								delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
								return ret;
							}							
						}
					  
					} else {
						if (verbosity >= 2)
							fprintf(stdout,	"# dropped region list covering only %i bases\n", num_read_map);
					}
				  
					current_seq = str;
					current_regions.clear();
					current_positions.clear();
				  
				  
					for (size_t i = 0; i < read.length(); i++)
						read_map[i] = regions[ori][chrN][nregion]->read_map[i];
				
					hits.topAlignments().update_top_alignment_indice();
					
				} else {
					// Regions are close enough that they may contain parts of one
					// spliced hit. They need thus be part of the same alignment
					// run. String this region onto the current sequence.
				  
					current_seq.append(str);
				  
					// merge read_maps
					for (size_t i = 0; i < read.length(); i++)
						read_map[i] = read_map[i] || regions[ori][chrN][nregion]->read_map[i];
				}
				current_regions.push_back(regions[ori][chrN][nregion]);
			  
				for (int p = 0; p < regions[ori][chrN][nregion]->end - regions[ori][chrN][nregion]->start; p++) 
				{
					current_positions.push_back(regions[ori][chrN][nregion]->start + p);
				}
			  
				/*fprintf(stdout, "read_map: ") ;
				  for (int i=0; i<_read.lenght(); i++)
				  if (read_map[i])
				  fprintf(stdout, "1") ;
				  else
				  fprintf(stdout, "0") ;
				  fprintf(stdout, "\n") ;*/
			}
		  
			int num_read_map = 0;
			for (size_t i = 0; i < read.length(); i++)
				if (read_map[i])
					num_read_map++;
		  
			if (verbosity >= 2)
				fprintf(stdout, "# region list covering %i bases\n", num_read_map);
			if (num_read_map >= _config.SPLICED_HIT_MIN_LENGTH_COMB)// && current_regions.size()>=2)
			{
				bool isunspliced ;
				//Recover long regions (starting points in sequences) for current regions to align
				std::vector<region_t*> corres_long_regions;
				corres_long_regions.clear();
				//fprintf(stdout,"chr length %i\n",genome->chromosome(chrN).length());
				recover_long_regions(read, corres_long_regions, long_regions[ori][chrN], current_regions);
				assert(corres_long_regions.size()>0); // at least one long region as support of alignment
				assert(corres_long_regions[0]->read_map!=NULL);

				//Take the first long region  to start alignment
				int hit_read_position = get_first_read_map(read, corres_long_regions[0]->read_map);
				int hit_len= corres_long_regions[0]->end-corres_long_regions[0]->start;
				if(ori==1){
					hit_read_position = read.length()-hit_len-hit_read_position+1;
				}
				assert (hit_read_position>=0 && hit_len >0);
				// fprintf(stdout,"read id %s curr len %i\n",read.id(), current_positions.size());
				// fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",hit_read_position, 
				//    corres_long_regions[0]->start, hit_len);
				// fprintf(stdout,	"# Number of current regions %i\n",(int)current_regions.size());					
				

				int transcription_direction = get_transcription_direction(_config.STRAND,ori);
						
				if (transcription_direction >=0)
				{
					int ret = perform_alignment_starter(result, hits, read_seq[ori], read_quality[ori],
														current_seq, current_regions, current_positions, chr, '+', ori,hit_read_position,
														corres_long_regions[0]->start, hit_len, non_consensus_search, num_alignments_reported,false,annotatedjunctions); 
					if (ret < 0)
					{
						result.delete_regions();
						delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
						return ret;
					}
			
					isunspliced = ret;
				}
				//if (!isunspliced) 
				if (transcription_direction <=0)
				{
					//fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",read.length()-(hit_read_position+hit_len),corres_long_regions[0]->end, hit_len);					  
					// fprintf(stdout,	"# Number of current regions %i\n",(int)current_regions.size());					  
					int ret = perform_alignment_starter(result, hits, read_seq[1 - ori],
														read_quality[1 - ori], current_seq,
														current_regions, current_positions, chr, '-', 1-ori,read.length()-(hit_read_position+hit_len),
														corres_long_regions[0]->end-1, hit_len, non_consensus_search, num_alignments_reported,false,annotatedjunctions);//end nucleotide in dna not included
					if (ret < 0)
					{
						result.delete_regions();
						delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
						return ret;
					}
				}
	  
				hits.topAlignments().update_top_alignment_indice();

			} else {
				if (verbosity >= 2)
					fprintf(stdout,	"# dropped region list covering only %i bases\n", num_read_map);
			}
		}
	}
	
	delete[] read_map;

	if (!_config.MAP_JUNCTIONS){
		result.delete_regions();
		delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
	}
	
	_stats.qpalma_region_align_time += clock() - start_time;
	
	if (verbosity >= 1 || ((clock()-_stats.qpalma_last_timing_report)/CLOCKS_PER_SEC>=10))
    {
		_stats.qpalma_last_timing_report = clock() ;
		_stats.qpalma_timing(((float) clock() - start_time) / CLOCKS_PER_SEC);
    }

	//fprintf(stderr, "num_alignments_reported=%i\n", ret) ;
	return num_alignments_reported ;

}


void *perform_alignment_wrapper(QPalma::perform_alignment_t *data)
{
	try
	{
		assert(data->qpalma!=NULL) ;
		data->ret = data->qpalma->perform_alignment(*data->result, *data->readMappings, data->read_string, data->read_quality, data->dna,
													data->current_regions, data->positions, *data->contig_idx,
													data->strand, data->ori, data->hit_read,
													data->hit_dna,data->hit_length, data->non_consensus_search, data->aln, data->remapping,*data->annotatedjunctions) ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr, "thread exception caught\n") ;
		data->ret = -1 ; // failed
		
	}

	return data ;
}


int QPalma::junctions_remapping(Hits &hits, Result &result, JunctionMap &junctionmap, int nb_spliced_alignments, JunctionMap &annotatedjunctions) const {

	std::vector<std::vector<region_t *> > *long_regions = result.long_regions;
	// if (nb_spliced_alignments > 0){
	// 	result.delete_regions();
	// 	delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
	// 	return 0;
	// }
	const int junction_tol = 10 ; // the hit may overlap by this length with the junction
	
	Read const &read(hits.getRead());

   	
	// Current regions correspond to regions around the junction
	std::vector<region_t *> current_regions;
	//Positions are real positions corresponding to current regions
	std::vector<int> current_positions;
	
	int num_alignments_reported = 0 ;
	
 	std::string read_seq[2];
	read_seq[0] = std::string(read.data(), read.length());
	read_seq[1] = reverse(complement(read_seq[0]));
	if (verbosity >= 3) {
		fprintf(stdout, "# read[0]: %s\n", read_seq[0].c_str());
		fprintf(stdout, "# read[1]: %s\n", read_seq[1].c_str());
	}
  
	std::string read_quality[2];
	read_quality[0] = std::string(read.quality(0), read.length());
	read_quality[1] = reverse(read_quality[0]);
	if (verbosity >= 3)
		fprintf(stdout, "# readqual[0]: %s\n", read_quality[0].c_str());
  

	for (int ori = 0; ori < 2; ori++){

		int transcription_direction = get_transcription_direction(_config.STRAND,ori);
						
		hits.topAlignments().init_top_alignment_indice();
	
		for (size_t chrN = 0; chrN < long_regions[ori].size(); chrN++) 
		{

			if (long_regions[ori][chrN].size() == 0)
				continue;

			current_regions.clear();
			current_positions.clear();
			
			Chromosome const &chr = genome->chromosome(chrN);			

			//Sort long regions by start positions
			region_t ** arr = NULL ;
			//size_t nbr_long_regions=long_regions.size();
			try 
			{
				arr = new region_t*[long_regions[ori][chrN].size()] ;			
			}
			catch (std::bad_alloc&)
			{
				fprintf(stderr, "[remapping_junctions] ERROR Could not allocate memory (_read.id() = %s)\n",
						read.id());
			}
			for (int i = 0; i < (int)long_regions[ori][chrN].size(); i++)
				arr[i] = long_regions[ori][chrN][i];  
			qsort(arr, long_regions[ori][chrN].size());

			for (size_t nregion = 0; nregion < long_regions[ori][chrN].size(); nregion++)
			{
				int rstart_=arr[nregion]->start;
				int rend_=arr[nregion]->end;
				
				if (verbosity>4)
					fprintf(stdout, "long region: %i-%i (%ld,%i)\n", rstart_, rend_, chrN, ori);

				// find a lower bound on the index with binary search
				junctionmap.lock() ;
				//std::deque<Junction>::iterator it = junctionmap.junctionlist[chrN].begin();
				std::deque<Junction>::iterator it = my_lower_bound(junctionmap.junctionlist[chrN].begin(), junctionmap.junctionlist[chrN].end(), rstart_-1000) ;
				junctionmap.unlock() ;
				
				while (it != junctionmap.junctionlist[chrN].end())
				{
					int rstart=rstart_ ;
					int rend=rend_ ;

					//Search for an overlapping with junction

					junctionmap.lock() ;

					//Intervals around junctions
					int int1_start = (*it).start - (read.length()+1);
					if (int1_start <0)
						int1_start=0;
					int int1_end = (*it).start-1;
					int int2_start = (*it).end + 1;					
					int int2_end = (*it).end + (read.length()+1);
					char strand = (*it).strand ;

					junctionmap.unlock() ;
					
					
					//Continue only if the strand with the splice junction is consistent with the transcription direction
					if ((transcription_direction==1 && strand=='-') || (transcription_direction==-1 && strand=='+') )		
						break;
					
					if ((unsigned int)int2_end >= chr.length())
						int2_end=chr.length()-1;



					if (rstart < int1_start)
						break; //no overlapping -> next long region
					
					//Overlapping junction
					if ( (rstart >= int1_start and rend <= int1_end+junction_tol+1) or (rstart >= int2_start-junction_tol and rend <= int2_end+1))
					{
						if (rstart >= int1_start and !(rend <= int1_end+1) and (rend <= int1_end+junction_tol+1)){
							
							rend=int1_end+1 ;

						}
						
						if (!(rstart >= int2_start) and rend <= int2_end+1 and (rstart >= int2_start-junction_tol)){
							
							rstart=int2_start ;
						}
						
					 	if (verbosity>4)
						{
							fprintf(stdout, "try to align\n") ;
							fprintf(stdout, "    junction: (%i-%i)-(%i-%i) (%ld,%c)\n",int1_start,int1_end,int2_start,int2_end,chrN,strand);
							fprintf(stdout, "    junction %c%c-%c%c\n",chr[(*it).start],chr[(*it).start+1],chr[(*it).end-1],chr[(*it).end]);
							fprintf(stdout, "long region: %i-%i (%ld,%i)\n", rstart, rend, chrN, ori);
						}

						//Create regions from junction
						region_t *new_region1 = NULL; 
						try {
							new_region1 = new region_t();
							new_region1->read_map = new bool[read.length()];
						} catch (std::bad_alloc&) 
						{
							fprintf(stderr, "[remapping_junctions] allocating memory for read_map failed\n");
							result.delete_regions();
							delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
							return -1;
						}
						new_region1->start = int1_start;
						new_region1->end = int1_end+1;
						new_region1->from_map = true ;
						for (size_t ii = 0; ii < read.length(); ii++)
							new_region1->read_map[ii] = false;
						
						region_t *new_region2 = NULL; 
						try {
							new_region2 = new region_t();
							new_region2->read_map = new bool[read.length()];
						} catch (std::bad_alloc&) 
						{
							fprintf(stderr, "[remapping_junctions] allocating memory for read_map failed\n");
							result.delete_regions();
							delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
							return -1;
						}
						new_region2->start = int2_start;
						new_region2->end = int2_end+1;
						new_region2->from_map = true ;
						for (size_t ii = 0; ii < read.length(); ii++)
							new_region2->read_map[ii] = false;
						

						current_regions.push_back(new_region1);
						current_regions.push_back(new_region2);
						
						//Corresponding positions and sequence

 						std::string current_seq;
						current_seq.assign("") ;
						for (int p = int1_start; p <= int1_end; p++) 
						{
							current_positions.push_back(p);
							current_seq.push_back(chr[p]) ;							
						}
						for (int p = int2_start; p <= int2_end; p++) 
						{
							current_positions.push_back(p);
							current_seq.push_back(chr[p]) ;							
						}


						//Take the first long region  to start alignment
						int hit_read_position = get_first_read_map(read, arr[nregion]->read_map);
						if (ori==0)
							hit_read_position += (rstart - rstart_);
						else
							hit_read_position += (rend_-rend);

						int hit_len= rend - rstart ; //arr[nregion]->end-arr[nregion]->start;


						if(ori==1){
							hit_read_position = read.length()-hit_len-hit_read_position +1;
						}

						assert (hit_read_position>=0 && hit_len >0);
						if (verbosity>3)
						{
							fprintf(stdout,"read id %s curr len %i\n",read.id(), (int)current_positions.size());
							fprintf(stdout,	"# Starting point for alignments: read %i, dna %i, len %i\n",hit_read_position, rstart /*arr[nregion]->start*/, hit_len);					  
							fprintf(stdout,	"# Number of current regions %i\n",(int)current_regions.size());
							
							fprintf(stdout,"DNA:%s\n",current_seq.c_str());
							fprintf(stdout,"READ:%s\n",read_seq[ori].c_str());
						}
						int ret;
						


						if (strand == '+')
						{
							ret = perform_alignment_starter(result, hits, read_seq[ori], read_quality[ori], 
															current_seq, current_regions, current_positions, 
															chr, '+', ori, hit_read_position, 
															rstart /*arr[nregion]->start*/, 
															hit_len, false, num_alignments_reported, true,
															annotatedjunctions);
						}
						else
						{
							ret = perform_alignment_starter(result, hits, read_seq[1 - ori], read_quality[1 - ori], 
															current_seq, current_regions, current_positions, 
															chr, '-', 1-ori, read.length()-(hit_read_position+hit_len),
															rend /*arr[nregion]->end*/ -1, hit_len, false, num_alignments_reported, true, annotatedjunctions);
						}
						
						if (ret < 0)
						{
							delete_long_regions(long_regions); 
							return ret;
						}
						
						num_alignments_reported += ret;						
						current_regions.clear();
						current_positions.clear();
				  						
						hits.topAlignments().update_top_alignment_indice();

					}//end alignment

					it++;
					
				}//loop through junctions
			}//loop through long regions for a chromosome and ori 
		}// loop chr
	}// loop ori
	
			
	result.delete_regions();
	delete_long_regions(long_regions); //Need to be deleted because of deep copies of region_t elements
	
	
	//fprintf(stdout, "num_alignments_reported=%i\n", num_alignments_reported) ;
	return num_alignments_reported ;
	
}


 

// TODO: dd remove relicts from multithreading
int QPalma::perform_alignment_starter(Result &result, 
									  Hits &readMappings, std::string read_string, std::string read_quality, std::string dna, std::vector<region_t *> current_regions, std::vector<int> positions, 
									  Chromosome const &contig_idx, char strand, int ori,int hit_read_position, int hit_dna_position, int hit_length, bool non_consensus_search, int& num_alignments_reported, bool remapping,JunctionMap &annotatedjunctions) const
{
	struct perform_alignment_t* data = NULL ;
	try
	{

		ALIGNMENT* non_consensus_alignment = NULL ;
		ALIGNMENT* consensus_alignment = NULL ;


		{
			data = new struct perform_alignment_t ;
			
			data->result = &result;
			data->readMappings = &readMappings;
			data->read_string=read_string ;
			data->read_quality=read_quality ;
			data->dna = dna ;
			data->current_regions=current_regions ;
			data->positions=positions ;
			data->contig_idx=&contig_idx;
			data->strand =strand ;
			data->ori = ori ;
			data->ret = -1000 ;
			data->hit_read = hit_read_position;
			data->hit_dna = hit_dna_position;
			data->hit_length = hit_length;
			data->qpalma = this ;
			data->joined=false ;
			data->aln=NULL ;
			data->non_consensus_search=false ;
			data->remapping=remapping;
			data->annotatedjunctions=&annotatedjunctions;
			

			perform_alignment_wrapper(data);

			consensus_alignment = data->aln ;
			data->aln = NULL ;
			data->joined = true;
			delete data ;
		}

		
		if (non_consensus_search) // && !consensus_alignment)
		{
			data = new struct perform_alignment_t ;

			data->result = &result;
			data->readMappings = &readMappings;
			data->read_string=read_string ;
			data->read_quality=read_quality ;
			data->dna = dna ;
			data->current_regions=current_regions ;
			data->positions=positions ;
			data->contig_idx=&contig_idx;
			data->strand =strand ;
			data->ori = ori ;
			data->ret = -1000 ;
			data->hit_read = hit_read_position;
			data->hit_dna = hit_dna_position;
			data->hit_length = hit_length;
			data->qpalma = this ;
			data->joined=false ;
			data->aln=NULL ;
			data->non_consensus_search=true ;
			data->remapping=remapping;
			data->annotatedjunctions=&annotatedjunctions;

			perform_alignment_wrapper(data);

			non_consensus_alignment=data->aln ;
			data->aln=NULL;

			delete data ;
		}

		
		
		if (consensus_alignment && (consensus_alignment->passed_filters || non_consensus_search))
		{
			readMappings.topAlignments().add_alignment_record(consensus_alignment, 1) ;
			num_alignments_reported++ ;
		}
		else
		{
			if (non_consensus_alignment && non_consensus_alignment->passed_filters)
			{
				readMappings.topAlignments().add_alignment_record(non_consensus_alignment, 1) ;
				num_alignments_reported++ ;
			}
			
		}
		
		
		return 0 ; // in the meaning of perform_alignment this corresponds to a spliced alignment
	}
	catch (std::bad_alloc&)
	{
	    fprintf(stdout, "WARNING: Thread creating failed\n") ;
	    delete data ;
	    
	    return -1 ;
	}
}


//TODO: remove me 08.06.2010
//int QPalma::rescue_alignment(Read const &read, std::string & read_anno, int ori, int &num_A, int &num_T, int &num)
//{
//	unsigned int read_pos = 0 ;
//	int genome_pos = 0 ;
//	int last_good_pos = 0 ;
//	int alignment_gaps = 0 ;
//	int alignment_mismatches = 0 ;
//
//  for (unsigned int i=0; i<read_anno.length(); i++)
//    {
//      //fprintf(stdout, "%c", read_anno[i]) ;
//      assert(read_anno[i]!=']') ;
//      if (read_anno[i]!='[')
//	{
//	  genome_pos++ ;
//	  read_pos++ ;
//	  num++ ;
//	  if (read_anno[i]=='A' || read_anno[i]=='a')
//	    num_A++ ;
//	  if (read_anno[i]=='T' || read_anno[i]=='t')
//	    num_T++ ;
//	}
//      else
//	{
//	  if (read_anno[i+1]=='-')
//	    {
//	      alignment_gaps++ ;
//	      if (ori==0)
//		read_pos++ ;
//	      else
//		genome_pos++ ;
//	      num++ ;
//	      if (read_anno[i+2]=='A' || read_anno[i+2]=='a')
//		num_A++ ;
//	      if (read_anno[i+2]=='T' || read_anno[i+2]=='t')
//		num_T++ ;
//	    }
//	  else if (read_anno[i+2]=='-')
//	    {
//	      if (ori==0)
//		genome_pos++ ;
//	      else
//		read_pos++ ;
//	      alignment_gaps++ ;
//	    }
//	  else
//	    {
//	      read_pos++ ;
//	      genome_pos++ ;
//	      alignment_mismatches++ ;
//	      num++ ;
//	      if (read_anno[i+2]=='A' || read_anno[i+2]=='a')
//		num_A++ ;
//	      if (read_anno[i+2]=='T' || read_anno[i+2]=='t')
//		num_T++ ;
//	    }
//	  i+=3 ;
//	}
//      if (!(read_pos<=read.length()))
//	fprintf(stderr, "ASSERT: %i, %i, %i, %i: %s\n", i, ori, read_pos, read.length(), read_anno.c_str()) ;
//
//      double frac=read_pos/read.length() ;
//      if (alignment_mismatches <= _config.NUM_MISMATCHES*frac && alignment_gaps <= _config.NUM_GAPS*frac && alignment_mismatches+alignment_gaps <= hits.*frac)
//	{
//	  last_good_pos = read_pos ;
//	  num_A=0 ;
//	  num_T=0 ;
//	  num=0 ;
//	}
//    }
//  return last_good_pos ;
//}


int QPalma::perform_alignment(Result &result, Hits &readMappings, std::string &read_string, std::string &read_quality, std::string &dna, std::vector<region_t *> &current_regions, std::vector<int> &positions, 
							  Chromosome const &contig_idx, char strand, int ori, int hit_read, int hit_dna, int hit_length, bool non_consensus_search, ALIGNMENT *& aln, bool remapping, JunctionMap &annotatedjunctions) const
// ori = read orientation
// strand = dna strand/orientation

//hit_read, hit_dna, hit_length: starting (real) positions and length of the hit for starting the alignment

{
	Read const &read(result._read);
	bool isunspliced = false; // is only true, if a valid alignment was reported and it is unspliced

	if (verbosity>=3)
	{
		for (size_t i=0; i<current_regions.size(); i++)
			fprintf(stdout, "region %i: %i - %i (%i)\n", (int)i, current_regions[i]->start, current_regions[i]->end, current_regions[i]->from_map) ;
	}
	
	if (verbosity>=2)
		fprintf(stdout, "average alignment length %lint\n", _stats.qpalma_total_dna_length/(_stats.qpalma_total_alignments+1)) ;
	
	// TODO
	//current_regions[0]->strand = strand ;
	//current_regions[0]->chromosome = contig_idx ;

	if (verbosity>=1)
		fprintf(stdout, "# [perform_alignment] performing alignment of read %s with %i regions and sequence of length %i, strand=%c\n", read.id(), (int)current_regions.size(), (int)dna.size(), strand) ;

	int nr_paths_p = 1;

	if (verbosity >= 3) {
		for (size_t i = 0; i < current_regions.size(); i++) {
			fprintf(stdout, "# region %i: %i - %i\n", (int)i, current_regions[i]->start,
					current_regions[i]->end);
		}
	}

	/* setup read sequence and quality information */
	int est_len_p = read_string.length();
	//char* est = (char *) malloc(sizeof(char) * (est_len_p + 1));
	char* est = NULL ;
	try
	{
		est = new char[est_len_p + 1] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr,	"[perform_alignment] Could not allocate memory (_read.id() = %s; est)\n", read.id());
		return -1;
	}
	strncpy(est, read_string.c_str(), est_len_p);
	est[est_len_p] = '\0';

	double* prb = NULL ;
	try
	{
		prb = new double[est_len_p] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr,	"[perform_alignment] Could not allocate memory (_read.id() = %s; prb)\n", read.id());
		delete[] est;
		return -1;
	}

	int offset_switched=0;
	for (size_t i = 0; i < (size_t) est_len_p; i++)
	{
		prb[i] = (read_quality[i] - alignment_parameters->quality_offset);
		if (prb[i]<-10 || prb[i]>70)
			fprintf(stderr, "prb[%i]=%f (offset=%i, %s, %s)\n", (int)i, prb[i], alignment_parameters->quality_offset, read_quality.c_str(), read.data()) ;

		if (_config.QPALMA_PRB_OFFSET_FIX)
		{
			if (alignment_parameters->quality_offset==33 && prb[i]>70)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 64) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=64 ;
				read.set_quality_offset(64) ;
				offset_switched=1;
				
			}
			if (alignment_parameters->quality_offset==64 && prb[i]<-10)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 33) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=33 ;
				read.set_quality_offset(33) ;
				offset_switched=1;
			}
		}
	}
	
	if (offset_switched==1){
		for (int i = 0; i < est_len_p; i++)
			prb[i] = (read_quality[i] - alignment_parameters->quality_offset);	
	}

	
	int num_qual_support = 0;
	if (alignment_parameters->num_qualityPlifs > 0)
		num_qual_support = alignment_parameters->qualityPlifs[0].len;
	//Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, true, verbosity);
	Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, !(_config.NO_QPALMA), verbosity);
	/* initialize acceptor and donor tables */
	int d_len = dna.length() ;
	double* donor = NULL ;
	try
	{
		donor = new double[d_len] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr, "[perform_alignment] Could not allocate memory (_read.id() = %s; donor=%i)\n", read.id(), d_len);
		return -1;
	}

	for (int i = 0; i < d_len; i++)
		donor[i] = -ALMOST_INFINITY;

	int a_len = dna.length();
	double* acceptor = NULL ;
	try
	{
		acceptor = new double[a_len] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr,	"[perform_alignment] Could not allocate memory (_read.id() = %s; acceptor)\n", read.id());
		delete[] donor;
		delete[] est ;
		delete[] prb;
		return -1;
	}
	for (int i = 0; i < a_len; i++)
		acceptor[i] = -ALMOST_INFINITY;

	const int num_intervals = current_regions.size();
	int interval_matrix[num_intervals * 2];
	size_t cum_length[num_intervals + 1];
	cum_length[0] = 0;
	for (int i = 0; i < num_intervals; i++) {
		interval_matrix[i * 2 + 0] = current_regions[i]->start;
		interval_matrix[i * 2 + 1] = current_regions[i]->end;
		cum_length[i + 1] = cum_length[i] + current_regions[i]->end - current_regions[i]->start;
	}

	
	assert(cum_length[num_intervals]==dna.length());

	

	if (!_config.NO_SPLICE_PREDICTIONS)
	{


		if (!remapping){
			
			if (_config.ACC_FILES.length()>0 && _config.DON_FILES.length()>0){
				
				/* query acceptor and donor scores */
				IntervalQuery iq;
				char basename[1000];
				const int num_scores = 1;
				const char *score_names[num_scores] = { "Conf_cum" };
	  

				// acc
				sprintf(basename, _config.ACC_FILES.c_str(), contig_idx.nr()+1, strand) ;
	    
				int acc_size = iq.query(basename, (char**) score_names, num_scores,
										interval_matrix, num_intervals);
				if (acc_size < 0) 
				{
					delete[] acceptor;
					delete[] donor ;
					delete[] est;
					delete[] prb;
					return -1;
				}
	    
				int* acc_pos = NULL ;
				int* acc_index = NULL ;
				double* acc_score = NULL ;
				try
				{
					acc_pos = new int[acc_size] ;
					acc_index = new int[acc_size] ;
					acc_score = new double[acc_size] ;
				}
				catch (std::bad_alloc&)
				{
					fprintf(stderr, "[perform_alignment] Could not allocate memory (_read.id() = %s)\n", read.id());
					delete[] donor;
					delete[] est;
					delete[] prb;
					delete[] acc_pos;
					delete[] acc_index;
					delete[] acc_score;
					return -1;
				}
	    
				iq.getResults(acc_pos, acc_index, acc_score);
				iq.cleanup();

				// don
				sprintf(basename, _config.DON_FILES.c_str(), contig_idx.nr()+1, strand) ;
	    
				int don_size = iq.query(basename, (char**) score_names, num_scores,
										interval_matrix, num_intervals);
				if (don_size < 0) {
					/* cleanup */
					delete[] acceptor;
					delete[] donor;
					delete[] est;
					delete[] prb;
					delete[] acc_pos;
					delete[] acc_index;
					delete[] acc_score;
					return -1;
				}
	  
				int* don_pos = NULL ;
				int* don_index = NULL ;
				double* don_score = NULL ;
				try
				{
					don_pos = new int[don_size] ;
					don_index = new int[don_size] ;
					don_score = new double[don_size] ;
				}
				catch (std::bad_alloc&)
				{
					fprintf(stderr, "[perform_alignment] Could not allocate memory (_read.id() = %s)\n", read.id());
					/* cleanup */
					delete[] donor ;
					delete[] est ;
					delete[] prb ;
					delete[] acc_pos ;
					delete[] acc_index ;
					delete[] acc_score ;
					delete[] don_pos ;
					delete[] don_index ;
					delete[] don_score ;
					return -1;
				}
				iq.getResults(don_pos, don_index, don_score);
				iq.cleanup();
		
		

// 	  /* write acceptor and donor predictions into tables */
// 	  if (strand=='+') // only necessary for + strand
// 	    for (int i=0; i<acc_size; i++)
// 	      {
// 		//fprintf(stdout,"%i ",acc_pos[i]);
// 		//if (strand=='+') // only necessary for + strand
// 		acc_pos[i]-=2 ; // 1-based -> 0-based  and shift by 1
// 		for (size_t j=0; j<current_regions.size(); j++)
// 		  if (acc_pos[i]>=current_regions[j]->start && acc_pos[i]<current_regions[j]->end){
// 		    acceptor[ acc_pos[i] - current_regions[j]->start + cum_length[j] ] = acc_score[i] ;
// 		  }
// 	      }
// 	  else
// 	    for (int i=0; i<acc_size; i++)
// 	      {
// 		//		fprintf(stdout,"%i ",acc_pos[i]);
// 		for (size_t j=0; j<current_regions.size(); j++)
// 		  if (acc_pos[i]>=current_regions[j]->start && acc_pos[i]<current_regions[j]->end)
// 		    acceptor[ acc_pos[i] - current_regions[j]->start + cum_length[j] ] = acc_score[i] ;
// 	      }
// 	  //	  fprintf(stdout,"\n");
		
// 	  for (int i = 0; i < don_size; i++) {
// 	    don_pos[i] -= 1; // 1-based -> 0-based
// 	    for (size_t j = 0; j < current_regions.size(); j++)
// 	      if (don_pos[i] >= current_regions[j]->start && don_pos[i]  < current_regions[j]->end)
// 		donor[don_pos[i] - current_regions[j]->start + cum_length[j]] = don_score[i];
// 	  }


				size_t j;
				size_t nbr_regions=current_regions.size();
				j=0;
	  
				//Splice sites and current regions are sorted by start position
				if (strand=='+'){
					for (int i=0; i<acc_size; i++){
						acc_pos[i]-=2 ; // 1-based -> 0-based  and shift by 1
						//fprintf(stdout,"(%i %i %i %i) ",acc_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						while (j<nbr_regions-1 && acc_pos[i]>=current_regions[j]->end){
							j++;
						}
	      
						//fprintf(stdout,"c(%i %i %i %i) ",acc_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (acc_pos[i]>=current_regions[j]->start && acc_pos[i]<current_regions[j]->end)
							acceptor[ acc_pos[i] - current_regions[j]->start + cum_length[j] ] = acc_score[i] ;  
					}
					j=0;
					for (int i = 0; i < don_size; i++) {
						don_pos[i] -= 1; // 1-based -> 0-based
						//fprintf(stdout,"b(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						while (j<nbr_regions-1 && don_pos[i] >= current_regions[j]->end)
							j++;
						//fprintf(stdout,"(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (don_pos[i]>=current_regions[j]->start && don_pos[i] < current_regions[j]->end){
							donor[don_pos[i] - current_regions[j]->start + cum_length[j]]= don_score[i];
						}
					}
				}
		
		
				
				//Negative strand: directly modify the relative coordinate of the splice site in reversed and complemented dna
				else{
					for (int i=0; i<acc_size; i++){
						while (j<nbr_regions-1 && acc_pos[i]>=current_regions[j]->end){
							j++;
						}

						//fprintf(stdout,"(%i %i %i %i) ",acc_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (acc_pos[i]>=current_regions[j]->start && acc_pos[i]<current_regions[j]->end){
							acceptor[a_len-1-(acc_pos[i] - current_regions[j]->start + cum_length[j]) ] = acc_score[i] ;
						}
					}
					//	  fprintf(stdout,"\n");
					j=0;
					for (int i = 0; i < don_size; i++) {
						don_pos[i] -= 1; // 1-based -> 0-based
						//fprintf(stdout,"b(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						while (j<nbr_regions-1 && don_pos[i] >= current_regions[j]->end)
							j++;
						//fprintf(stdout,"(%i %i %i %i) ",don_pos[i],current_regions[j]->start,current_regions[j]->end,j);
						if (don_pos[i]>=current_regions[j]->start && don_pos[i] < current_regions[j]->end){
							donor[d_len-1 -(don_pos[i] - current_regions[j]->start + cum_length[j])]= don_score[i];
						}
					}
				}

				if (verbosity>=3)
					fprintf(stdout, "got %i acceptors and %i donors\n", acc_size, don_size) ;

				/* cleanup */
				delete[] acc_pos ;
				delete[] acc_index ;
				delete[] acc_score ;
				delete[] don_pos ;
				delete[] don_index ;
				delete[] don_score ;
			}
			
			//Use splice sites from the annotation
			if (_config.SCORE_ANNOTATED_SPLICE_SITES){

				//fprintf(stdout,"Use annotated splice sites\n");
				std::vector<int> acc_pos_vec;
				std::vector<int> don_pos_vec;

				//Get splice sites for different intervals				
				for (int i = 0; i < num_intervals; i++) {
					//fprintf(stdout,"Get splice sites for region %i\n",i);
					get_annotated_splice_sites(acc_pos_vec,don_pos_vec,annotatedjunctions,current_regions[i]->start,current_regions[i]->end,contig_idx.nr(),strand);
					//fprintf(stdout,"Get splice sites for region %i ---- DONE\n",i);
					//Transform coordinates
					if (strand=='+'){
						for (int k=0; k<(int)acc_pos_vec.size(); k++){
							//fprintf(stdout,"(%i %i %i %i)\n",acc_pos_vec[k],current_regions[i]->start,current_regions[i]->end,i);
							acceptor[ acc_pos_vec[k] - current_regions[i]->start + cum_length[i]] = 1.0 ;  
							//fprintf(stdout,"acceptor=%c|%c%c|%c\n",dna[acc_pos_vec[k] - current_regions[i]->start + cum_length[i]-2],dna[acc_pos_vec[k] - current_regions[i]->start + cum_length[i]-1],dna[acc_pos_vec[k] - current_regions[i]->start + cum_length[i]],dna[acc_pos_vec[k] - current_regions[i]->start + cum_length[i]+1]);
							
						}
						for (int k=0; k<(int)don_pos_vec.size(); k++) {
							donor[don_pos_vec[k] - current_regions[i]->start + cum_length[i]]= 1.0;
							//fprintf(stdout,"donor=%c|%c%c|%c\n",dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]-1],dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]],dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]+1],dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]+2]);
						}
					}

					else{
						
						for (int k=0; k<(int)acc_pos_vec.size(); k++){
							acceptor[a_len-1-(acc_pos_vec[k] - current_regions[i]->start + cum_length[i])] = 1.0 ;
						}
						for (int k=0; k<(int)don_pos_vec.size(); k++) {
							donor[d_len-1 -(don_pos_vec[k] - current_regions[i]->start + cum_length[i])]= 1.0;
							//fprintf(stdout,"donor=%c|%c%c|%c\n",dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]-2],dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]-1],dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]],dna[don_pos_vec[k] - current_regions[i]->start + cum_length[i]+1]);
						}
					}
					
					acc_pos_vec.clear();
					don_pos_vec.clear();
					
				}
			}
			

	

		}
	}
	
	if (strand=='-')
	{
		dna = reverse(complement(dna)) ;
		//reverse(donor, d_len) ;
		//reverse(acceptor, a_len) ;
		//for (int i=0; i<positions.size(); i++)
		//fprintf(stdout, "%i:%i ", i, positions[i]) ;
		//fprintf(stdout, "\n") ;
		positions = reverse(positions);
	}
	/*if (verbosity>=2)
	  {
	  fprintf(stdout, "positions:\n") ;
	  for (int i=0; i<positions.size(); i++)
	  fprintf(stdout, "%i:%i ", i, positions[i]) ;
	  fprintf(stdout, "\n") ;
	  }*/
	
	if (verbosity >= 3)
		fprintf(stdout, "# dna: %s\n", dna.c_str());
	if (verbosity >= 3)
		fprintf(stdout, "# read: %s\n", read_string.c_str());
	
	const float NON_CONSENSUS_SCORE = -123456;

	/* check whether we have scores for all donor and acceptor positions (first 10% of reads)*/
	if (!remapping){
		
		int match = 0, num = 0;
		for (int i = 2; i < d_len - 2; i++, num++)
		{
			bool is_ss = false ;
			for (unsigned int j=0; j<_config.DON_CONSENSUS.size(); j++)
				if (dna[i] == _config.DON_CONSENSUS[j][0] && dna[i+1] == _config.DON_CONSENSUS[j][1])
				{
					is_ss = true ;
					break ;
				}
			if (non_consensus_search && !is_ss) 
			{
				if (dna[i]!='N' && dna[i+1]!='N')
					donor[i] = NON_CONSENSUS_SCORE ;
				match += 1 ;//(donor[i] > -ALMOST_INFINITY);
			} 
			else
				if (is_ss) 
				{
					// if no splice predictions, put score 0 for allowed donor splice consensus (e.g. 'GT/C')
					if (_config.NO_SPLICE_PREDICTIONS)
						donor[i] = 0.0 ;
					if (_config.SCORE_ANNOTATED_SPLICE_SITES && donor[i]==-ALMOST_INFINITY)
						donor[i] = 0.0 ;
					match += (donor[i] > -ALMOST_INFINITY || _config.SCORE_ANNOTATED_SPLICE_SITES);
				} 
				else 
				{
					match += _config.SCORE_ANNOTATED_SPLICE_SITES || (donor[i] <= -ALMOST_INFINITY)  || dna[i]=='N' || dna[i+1]=='N' ;
					donor[i] = -ALMOST_INFINITY;
				}
		}
	
		if (match < num * 0.9)
			fprintf(stderr,
					"Warning: donor predictions do not match genome positions (match=%i  num=%i  strand=%c  ori=%c  chr=%s  start=%i  end=%i)\n",
					match, num, strand, ori == 0 ? '+' : '-', contig_idx.desc(),
					current_regions[0]->start, current_regions[current_regions.size() - 1]->end);
		//assert(match>=num*0.9) ; // otherwise positions will are shifted somehow
	

		match = 0; num = 0;
		for (int i = 2; i < a_len - 2; i++, num++)
		{
			bool is_ss = false ;
			for (unsigned int j=0; j < _config.ACC_CONSENSUS.size(); j++)
				if (i>0 && dna[i-1] == _config.ACC_CONSENSUS[j][0] && dna[i] == _config.ACC_CONSENSUS[j][1])
				{
					is_ss = true ;
					break ;
				}

			if (non_consensus_search && !is_ss)
			{
				if (dna[i]!='N' && dna[i-1]!='N')
					acceptor[i] = NON_CONSENSUS_SCORE ;
				match += 1;
			} 
			else
				if (is_ss)
				{
					if (_config.NO_SPLICE_PREDICTIONS)
						acceptor[i] = 0.0 ;
					if (_config.SCORE_ANNOTATED_SPLICE_SITES && acceptor[i]==-ALMOST_INFINITY)
						acceptor[i] = 0.0 ;
					match += (acceptor[i] > -ALMOST_INFINITY || _config.SCORE_ANNOTATED_SPLICE_SITES);
				} 
				else 
				{
					//if (acceptor[i]>-ALMOST_INFINITY && dna[i]!='N' && dna[i-1]!='N')
					//	fprintf(stdout, "acc over %i %f  %c%c\n", i, acceptor[i], dna[i-1], dna[i]) ;
					match += _config.SCORE_ANNOTATED_SPLICE_SITES || (acceptor[i] <= -ALMOST_INFINITY) || dna[i]=='N' || dna[i-1]=='N' ;
					acceptor[i] = -ALMOST_INFINITY;
				}
		}

		if (match<num*0.9)

			fprintf(stderr, "Warning: acceptor predictions do not match genome positions (match=%i  num=%i  strand=%c  ori=%c  chr=%s  start=%i  end=%i, %i)\n", 
					match, num, strand, ori==0 ? '+' : '-' , contig_idx.desc(), current_regions[0]->start, current_regions[current_regions.size()-1]->end, non_consensus_search) ;


		/* apply donor and acceptor plifs */
		for (int i = 0; i < d_len; i++)
			if (donor[i] > -ALMOST_INFINITY)
			{
				if (donor[i]==NON_CONSENSUS_SCORE)
				{
					donor[i]=0.01 ; // low donor probability  -> correspondingly low score
					if (not _config.NO_QPALMA)
						donor[i] = lookup_penalty(&alignment_parameters->d, 0, &donor[i]);
				}
				else
				{
					if (not _config.NO_QPALMA)
						donor[i] = lookup_penalty(&alignment_parameters->d, 0, &donor[i]);
				}			
			}
		for (int i = 0; i < a_len; i++)
			if (acceptor[i] > -ALMOST_INFINITY)
			{
				if (acceptor[i]==NON_CONSENSUS_SCORE)
				{
					acceptor[i]=0.01 ; // low acceptor probability  -> correspondingly low score
					if (not _config.NO_QPALMA)
						acceptor[i] = lookup_penalty(&alignment_parameters->a, 0, &acceptor[i]);
				}
				else
				{
					if (not _config.NO_QPALMA)
						acceptor[i] = lookup_penalty(&alignment_parameters->a, 0, &acceptor[i]);
				}
			}
	}
	
	bool remove_duplicate_scores = false;
//	bool print_matrix = false;

	clock_t start_time = clock();

	// Convert real dna position of the hit into relative dna position 
	int hit_dna_converted = convert_dna_position(hit_dna, cum_length, current_regions);
	if (strand=='-')
		hit_dna_converted=(int)dna.length()-1-hit_dna_converted;

	assert (hit_dna_converted >= 0);
	if (non_consensus_search)
	{
		alignment.myalign_fast(strand, contig_idx, positions, nr_paths_p, (char*) dna.c_str(), (int) dna.length(), est,
							   est_len_p, prb, alignment_parameters->h,
							   alignment_parameters->matchmatrix,
							   alignment_parameters->matchmatrix_dim[0]
							   * alignment_parameters->matchmatrix_dim[1], donor, d_len,
							   acceptor, a_len, alignment_parameters->qualityPlifs,
							   remove_duplicate_scores, hit_read, hit_dna_converted, hit_length, _config.SPLICED_MAX_INTRONS,
							   _config.NUM_GAPS, _config.NUM_MISMATCHES, readMappings.get_num_edit_ops(), 
							   MIN_NUM_MATCHES+ _config.MIN_NUM_MATCHES_PEN, remapping);
	}
	else
		alignment.myalign_fast(strand, contig_idx, positions, nr_paths_p, (char*) dna.c_str(), (int) dna.length(), est,
							   est_len_p, prb, alignment_parameters->h,
							   alignment_parameters->matchmatrix,
							   alignment_parameters->matchmatrix_dim[0]
							   * alignment_parameters->matchmatrix_dim[1], donor, d_len,
							   acceptor, a_len, alignment_parameters->qualityPlifs,
							   remove_duplicate_scores,hit_read,hit_dna_converted,hit_length,_config.SPLICED_MAX_INTRONS,
							   _config.NUM_GAPS, _config.NUM_MISMATCHES, readMappings.get_num_edit_ops(), 
							   MIN_NUM_MATCHES, remapping);

	static pthread_mutex_t clock_mutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_lock( &clock_mutex) ;
	_stats.qpalma_align_time += clock() - start_time;
	pthread_mutex_unlock( &clock_mutex) ;
	

	int s_align[dna.length()];
	int e_align[est_len_p];
	int mmatrix_p[alignment_parameters->matchmatrix_dim[0]
				  * alignment_parameters->matchmatrix_dim[1]];
	double alignscore;
	double *qScores = NULL;


	alignment.getAlignmentResults(s_align, e_align, mmatrix_p, &alignscore, qScores);
	int result_length = alignment.getResultLength();
	std::vector<int> exons;
	exons.clear();

	int exon_start = -1;
	int min_exon_len = read.length();
	int max_intron_len = 0;
	int min_intron_len = 1000000;

	bool alignment_valid=true ;

	int alignment_matches = 0;
	int alignment_gaps = 0;
	int alignment_mismatches = 0 ;
	int alignment_qual_mismatches = 0 ;
	std::string read_anno = std::string("");

	bool rescued_alignment=false;

	//Alignment not found with less than _config.NUM_EDIT_OPS or _config.NUM_GAPS or _config.NUM_MISMATCHES
	if (result_length<est_len_p)
	{
	    alignment_valid=false;
	    //fprintf(stdout, "result_length=%i, est_len_p=%i\n", result_length, est_len_p) ;
	    if (verbosity>=2)
			fprintf(stdout, "No alignment found\n") ;
	}
	else
	{  	
		
		bool intron_location=false;
		
	    for (size_t i = 0; i < dna.length(); i++) 
		{
			if (exon_start == -1 && s_align[i] == 0 && i<dna.length()-1) 
			{
				exon_start = i;
				continue;
			}
			if (exon_start!=-1 && i>0)
			{
				alignment_valid = remapping || (alignment_valid && ((positions[i-1]+1 == positions[i]) || (positions[i-1] == positions[i]+1))) ;
				intron_location=remapping && s_align[i]==0 && 
					((strand=='+' && positions[i-1]+1 != positions[i]) || 
					 (strand=='-' && positions[i-1] != positions[i]+1)) ;
				//fprintf(stdout, "pos[%i]=%i   pos[%i]=%i  s_align=%i strand=%c\n", (int)i-1, (int)positions[i-1], (int)i, (int)positions[i], s_align[i],strand) ;
				if (verbosity>=4)
					fprintf(stdout, "pos[%i]=%i   pos[%i]=%i  valid=%i\n", (int)i-1, (int)positions[i-1], (int)i, (int)positions[i], alignment_valid) ;
			}
			
			if (exon_start!=-1 && (s_align[i]!=0 || i==dna.length()-1 || intron_location))
			{
				if (exons.size()>0)
				{
					int intron_len = positions[exon_start]-exons[exons.size()-1] ;
					if (strand=='-')
						intron_len*=-1 ;
					if (intron_len>max_intron_len)
						max_intron_len=intron_len ;
					if (intron_len<min_intron_len)
						min_intron_len=intron_len ;
				}
				exons.push_back(positions[exon_start]) ;
				
				int exon_len;
				if (s_align[i]==0 && !intron_location){
					exons.push_back(positions[i]) ;
					//fprintf(stdout,"Exon: %i-%i\n",exon_start,i);
					exon_len=positions[i] - positions[exon_start] ;		
				}				
				else{					
					exons.push_back(positions[i-1]) ;
					//fprintf(stdout,"Exon: %i-%i\n",exon_start,i-1);
					exon_len=positions[i-1] - positions[exon_start] ;		
				}
								
				
				if (strand=='-')
					exon_len*=-1 ;
				exon_len++ ;
				assert(exon_len>=1) ;
				if (exon_len<min_exon_len)
					min_exon_len = exon_len ;
				if (intron_location)
					exon_start = i ;
				else
					exon_start = -1 ;
				
				intron_location=false;
				
				continue ;
			}
		}
	    assert(exon_start==-1) ;
	    assert(exons.size()>0) ;
	    if (verbosity>=2)
			fprintf(stdout, "alignment valid=%i (0 means that exons went over block boundaries)\n", alignment_valid) ;
	    
	   
	    if (strand=='-')
			exons=reverse(exons) ;
	    for (size_t i=0; i<exons.size(); i+=2)
			exons[i+1]++ ;
	    
	    if (verbosity >= 2) {
			fprintf(stdout, "# exons:\n");
			for (size_t i = 0; i < exons.size(); i += 2)
				fprintf(stdout, "# %i. %i - %i\n", (int)i / 2, exons[i], exons[i + 1]);
	    }

	    /*fprintf(stdout, "DNA: %d ", s_align[0]);
	      for (int i = 0; i < dna.length(); i++)
	      fprintf(stdout, "%i ", s_align[i]);
	      fprintf(stdout, "\n");
	      
	      fprintf(stdout, "EST: ");
	      for (int i = 0; i < est_len_p; i++)
	      fprintf(stdout, "%i ", e_align[i]);
	      fprintf(stdout, "\n");*/
	    
	    int dna_align[result_length];
	    int est_align[result_length];
	    char dna_align_str[result_length + 1];
	    char est_align_str[result_length + 1];
	    
	    alignment.getAlignmentArrays(dna_align, est_align);
	    
	    
	    //	int alignment_matches = 0;
	    //int alignment_gaps = 0;
	    //int alignment_mismatches = 0 ;
	    
	    //std::string read_anno = std::string("");
	    char map[8] = "-ACGTN*";
		int read_pos=0 ;
	    for (int i = 0; i < result_length; i++) 
		{
			assert(dna_align[i]>=0 && dna_align[i]<=6);
			dna_align_str[i] = map[dna_align[i]];
			
			assert(est_align[i]>=0 && est_align[i]<=6);
			est_align_str[i] = map[est_align[i]];
			
			if (est_align[i]!=0 && est_align[i]!=6 && dna_align[i]!=0)
			{
				if (est_align[i]==dna_align[i])
				{
					read_anno.push_back(map[est_align[i]]) ;
					if (map[est_align[i]]!='N')
						assert(map[est_align[i]]==read_string[read_pos]) ;
					read_pos++ ;
				}
				else
				{
					read_anno.push_back('[') ;
					if (strand=='+')
					{
						read_anno.push_back(map[dna_align[i]]) ;
						read_anno.push_back(map[est_align[i]]) ;
					} else
					{
						read_anno.push_back(map[est_align[i]]) ;
						read_anno.push_back(map[dna_align[i]]) ;
					}
					read_anno.push_back(']');
					if (map[dna_align[i]]!='N' && map[est_align[i]]!='N')
					{
						alignment_mismatches++ ;
						alignment_qual_mismatches += read_quality[read_pos]-read.get_quality_offset() ; 
						assert(map[est_align[i]]==read_string[read_pos]) ;
					}
					read_pos++ ;
				}
				alignment_matches += (est_align[i]==dna_align[i]) ;
			}
			else if ( est_align[i]==0 )
			{
				if (dna_align[i]!=0)
				{
					read_anno.push_back('[') ;
					if (strand=='+')
					{
						read_anno.push_back(map[dna_align[i]]) ;
						read_anno.push_back('-') ;
					}
					else
					{
						read_anno.push_back('-') ;
						read_anno.push_back(map[dna_align[i]]) ;
					}
					read_anno.push_back(']');
				}
				alignment_gaps++ ;
			}
			else if ( dna_align[i]==0)
			{
				read_anno.push_back('[') ;
				if (strand=='+')
				{
					read_anno.push_back('-') ;
					read_anno.push_back(map[est_align[i]]) ;
				}
				else
				{
					read_anno.push_back(map[est_align[i]]) ;
					read_anno.push_back('-') ;
				}
				read_anno.push_back(']');
				alignment_gaps++;
				if (map[est_align[i]]!='N')
					alignment_qual_mismatches += read_quality[read_pos]-read.get_quality_offset() ; 
				if (map[est_align[i]]!='N')
					assert(map[est_align[i]]==read_string[read_pos]) ;
				read_pos++ ;
			}
	    }
		assert(read_pos==(int)read.length()) ;
	    dna_align_str[result_length] = 0;
	    est_align_str[result_length] = 0;
	    
	    /*if (exons.size() == 2 && strand=='+' && ori==0)
	      {
	      fprintf(stderr, "unspliced alignscore=%2.4f\n", alignscore) ;
	      
	      float score2=score_unspliced(read_anno.c_str()) ;
	      //float score2=alignment.scoreUnsplicedAlignment(read_anno.c_str(), prb, _read.lenght(), alignment_parameters->qualityPlifs, alignment_parameters->matchmatrix) ;
	      fprintf(stderr, "recomputed alignscore=%2.4f\n", score2) ;
	      }*/


		if (strand=='-')
			read_anno=reverse(complement(read_anno)) ;

/*
  bool rescued_alignment = false ;
  int rescue_start = 0 ;
  int rescue_end = _read.lenght() ;
  if (alignment_valid && (alignment_mismatches > _config.NUM_MISMATCHES || alignment_gaps > _config.NUM_GAPS))
  {
  // reverse order of quality 
  char qual[500] ;
	    
  if (ori==0)
  for (int i=0; i<((int)_read.lenght()); i++)
  qual[i]=_read.quality()[0][i] ;
  else
  for (int i=0; i<((int)_read.lenght()); i++)
  qual[i]=_read.quality()[0][((int)_read.lenght())-i-1] ;
  qual[((int)_read.lenght())]=0 ;
	    
  int SL1_pos = -1 ;
  const int max_num_mismatch=2 ;
  const char *SL1="ggtttaattacccaagtttgag" ;
  for (int i=0; i<_read.lenght()-strlen(SL1); i++)
  {
  int num_mismatch = 0 ;
  for (int j=0; j<strlen(SL1); j++)
  {
  if (toupper(READ[i+j])!=toupper(SL1[j]))
  {
  num_mismatch++ ;
  if (num_mismatch>max_num_mismatch)
  break ;
  }
  }
  if (num_mismatch<=max_num_mismatch)
  {
  SL1_pos = i ;
  fprintf(stdout, "found SL1 at pos %i", SL1_pos) ;
  exit(-1) ;
  break ;
  }
  }
		
  int num_A=0, num_T=0, num_=0 ;
  if (ori==0)
  rescue_end = rescue_alignment(read_anno, ori, num_A, num_T, num_) ;
  else
  {
  std::string rev=reverse(complement(read_anno)) ;
  rescue_start = _read.lenght()-rescue_alignment(rev, ori, num_A, num_T, num_) ;
  }
  if (rescue_end-rescue_start>25)// && num_>0 && (num_A/num_>=0.5 || num_T/num_>=0.5) )//_read.lenght()/2)
  {
  //fprintf(stdout, "rescued poor alignment ori=%i, strand=%c, num_A=%i, num_T=%i: %s\t%s\n", ori, strand, num_A, num_T, read_anno.c_str(), qual) ;
  //fprintf(stdout, "rescued region: %i - %i \n", rescue_start, rescue_end) ;
  rescued_alignment = true ;
  }
  }
*/
		//	bool rescued_alignment = false ;

		if (verbosity >= 3) {
			fprintf(stdout, "# DNA : %s\n", dna_align_str);
			fprintf(stdout, "# Read: %s\n", est_align_str);
		}

		if (verbosity >= 2)
			fprintf(stdout,
					"# alignment: score=%1.3f  matches=%i  gaps=%i  anno=%s\n",
					alignscore, alignment_matches, alignment_gaps,
					read_anno.c_str());

		if (verbosity >= 1)
			fprintf(stdout, "# alignment with %i exons (%i, %i, %i, %i)\n", (int)exons.size()/2, 
					(int)alignment_mismatches, (int)alignment_gaps, (int)alignment_mismatches+alignment_gaps, (int)alignment_qual_mismatches) ;
	  
	}
	//if (alignment_matches >= read_string.length() - _config.NUM_EDIT_OPS
	//		&& exons.size() >= 4) // it must be spliced and not have too many mismatches

	bool alignment_passed_filters= (max_intron_len<_config.SPLICED_LONGEST_INTRON_LENGTH) && (min_intron_len>=_config.SPLICED_SHORTEST_INTRON_LENGTH) && 
		(rescued_alignment || (alignment_mismatches <= _config.NUM_MISMATCHES && alignment_gaps <= _config.NUM_GAPS && alignment_mismatches+alignment_gaps <= _config.NUM_EDIT_OPS))
		&&(exons.size()==2 ||( (exons.size() >= 4 || rescued_alignment) && ((int)exons.size() <= (_config.SPLICED_MAX_INTRONS+1)*2) && (min_exon_len >= _config.SPLICED_MIN_SEGMENT_LENGTH || remapping))) ;

	bool non_consensus_alignment=false ;

	if (alignment_valid) // it must be spliced and not have too many mismatches
	{
		aln = NULL;
		try 
		{
			aln = new ALIGNMENT();
		} 
		catch (std::bad_alloc&) 
		{
			fprintf(stderr,	"[capture_hits] allocating memory for aligment failed\n");
			result.delete_regions();
			return -1;
		}

		for (unsigned int ne=0; ne<(exons.size()/2)-1; ne++)
		{
			bool non_consensus_intron = false ;
			
		    int istart=exons[ne*2+1] ;
		    int istop =exons[(ne+1)*2] ;
		    char buf[1000] ;
		    if ((ori==0 && strand=='+') || (ori==1 && strand=='+'))
			{
				sprintf(buf, "%s:%i:%i:%c%c|%c%c%c%c..%c%c%c%c|%c%c", contig_idx.desc(), istart, istop-1, 
						contig_idx[istart-2], contig_idx[istart-1], contig_idx[istart], contig_idx[istart+1],  contig_idx[istart+2], contig_idx[istart+3], 
						contig_idx[istop-4], contig_idx[istop-3], contig_idx[istop-2], contig_idx[istop-1], contig_idx[istop], contig_idx[istop+1]) ;
				//fprintf(stdout, "intron+: %s\n", buf) ;
				
				bool is_consensus_don = false ;
				for (unsigned int j=0; j < _config.DON_CONSENSUS.size(); j++)
					if (contig_idx[istart] == _config.DON_CONSENSUS[j][0] && contig_idx[istart+1] == _config.DON_CONSENSUS[j][1])
					{
						is_consensus_don = true ;
						break ;
					}
				bool is_consensus_acc = false ;
				for (unsigned int j=0; j < _config.ACC_CONSENSUS.size(); j++)
					if (contig_idx[istop-2] == _config.ACC_CONSENSUS[j][0] && contig_idx[istop-1] == _config.ACC_CONSENSUS[j][1])
					{
						is_consensus_acc = true ;
						break ;
					}
				if (!is_consensus_don || !is_consensus_acc)
				{
					non_consensus_intron=true ;
					non_consensus_alignment=true ;
				}
				if (!non_consensus_search && !remapping)
					assert(!non_consensus_intron) ;
			}
		    else
			{
				sprintf(buf, "%s:%i:%i:%c%c|%c%c%c%c..%c%c%c%c|%c%c", contig_idx.desc(), istop-1, istart, 
						complement(contig_idx[istop+1]), complement(contig_idx[istop]), complement(contig_idx[istop-1]), complement(contig_idx[istop-2]), complement(contig_idx[istop-3]), complement(contig_idx[istop-4]), 
						complement(contig_idx[istart+3]), complement(contig_idx[istart+2]), complement(contig_idx[istart+1]), complement(contig_idx[istart-0]), complement(contig_idx[istart-1]), complement(contig_idx[istart-2])) ;
				//fprintf(stdout, "intron-: %s\n", buf) ;

				bool is_consensus_don = false ;
				for (unsigned int j=0; j < _config.DON_CONSENSUS.size(); j++)
					if (complement(contig_idx[istop-1]) == _config.DON_CONSENSUS[j][0] && complement(contig_idx[istop-2]) == _config.DON_CONSENSUS[j][1])
					{
						is_consensus_don = true ;
						break ;
					}
				bool is_consensus_acc = false ;
				for (unsigned int j=0; j < _config.ACC_CONSENSUS.size(); j++)
					if (complement(contig_idx[istart+1]) == _config.ACC_CONSENSUS[j][0] && complement(contig_idx[istart-0]) == _config.ACC_CONSENSUS[j][1]) 
					{
						is_consensus_acc = true ;
						break ;
					}
				if (!is_consensus_don || !is_consensus_acc)
				{
					non_consensus_intron=true ;
					non_consensus_alignment=true ;
				}
				if (!non_consensus_search && !remapping)
					assert(!non_consensus_intron) ;

			}
		    aln->non_consensus_intron.push_back(non_consensus_intron) ;
		    aln->intron_consensus.push_back(strdup(buf)) ;
		}

		if (remapping && exons.size()>=4){
			REAL score = 1.0 ; // high probability for donor and acceptor -> large ss score
			alignscore += lookup_penalty(&alignment_parameters->d, 0, &score) ; //alignment_parameters->d.penalties[alignment_parameters->d.len-1] ;
			alignscore += lookup_penalty(&alignment_parameters->a, 0, &score) ; //alignment_parameters->a.penalties[alignment_parameters->a.len-1] ;
		}
		
		aln->qpalma_score = alignscore;
		aln->num_matches = alignment_matches;
		aln->num_gaps = alignment_gaps;
        aln->num_mismatches = alignment_mismatches ;
        aln->qual_mismatches = alignment_qual_mismatches ;
		strcpy(aln->read_anno, read_anno.c_str());
		aln->exons = exons;
		aln->chromosome = &contig_idx;
		aln->strand = strand;
		//aln->rescued=rescued_alignment ;
		//aln->rescue_start = rescue_start ;
		//aln->rescue_end = rescue_end ;
		aln->passed_filters=alignment_passed_filters ;
		aln->non_consensus_alignment = non_consensus_alignment ; 
		aln->remapped = remapping ;

		aln->from_gm = 3;

		if (ori == 0)
			aln->orientation = '+';
		else
			aln->orientation = '-' ;
		aln->min_exon_len = min_exon_len ;
		aln->max_intron_len = max_intron_len ;
		aln->spliced = (exons.size()>=4) ;
		if (!aln->spliced)
			isunspliced = true;

		strcpy(aln->read_id, read.id());
		
		if (verbosity >= 2)
		{
			fprintf(stdout,
					"# alignment (%i/%i) with %i exons found for %s (score=%1.3f  matches=%i  gaps=%i strand=%c): %s\n",
					non_consensus_search,non_consensus_alignment, (int) exons.size() / 2, read.id(), alignscore,
					alignment_matches, alignment_gaps, strand, read_anno.c_str());
			for (size_t i = 0; i < exons.size(); i += 2)
				fprintf(stdout, "# exon %i: %i - %i\n", (int)i / 2, exons[i], exons[i+ 1]);
		}
	} 
	else
		if (verbosity>=1)
		{
			fprintf(stdout, "# dropped alignment with %i exons (%i, %i, %i)\n", (int)exons.size()/2, 
					(int)alignment_mismatches, (int)alignment_gaps, (int)alignment_mismatches+alignment_gaps) ;
			fprintf(stdout, "alignment_valid=%i\n", alignment_valid) ;
			fprintf(stdout, "max_intron_len<_config.SPLICED_LONGEST_INTRON_LENGTH=%i\n", max_intron_len<_config.SPLICED_LONGEST_INTRON_LENGTH) ;
			fprintf(stdout, "alignment_mismatches(%i) <= _config.NUM_MISMATCHES(%i)=%i\n", alignment_mismatches, _config.NUM_MISMATCHES, alignment_mismatches <= _config.NUM_MISMATCHES) ;
			fprintf(stdout, "alignment_gaps <= _config.NUM_GAPS=%i\n", alignment_gaps <= _config.NUM_GAPS) ;
			fprintf(stdout, "alignment_mismatches(%i)+alignment_gaps(%i) <= readMappings.get_num_edit_ops()=%i\n",alignment_mismatches, alignment_gaps, alignment_mismatches+alignment_gaps <= readMappings.get_num_edit_ops()) ;
			fprintf(stdout, "exons.size()=%i\n", (int)exons.size()) ;
			fprintf(stdout, "((int)exons.size() <= (_config.SPLICED_MAX_INTRONS+1)*2)=%i\n", ((int)exons.size() <= (_config.SPLICED_MAX_INTRONS+1)*2)) ;
			fprintf(stdout, "min_exon_len >= _config.SPLICED_MIN_SEGMENT_LENGTH=%i\n", min_exon_len >= _config.SPLICED_MIN_SEGMENT_LENGTH) ;
		}

	delete[] acceptor ;
	delete[] donor ;
	delete[] est ;
	delete[] prb ;

	return (int) isunspliced;
}

double QPalma::score_unspliced(Read const &read, const char * read_anno, const char strand, const char ori) const
{
	if (alignment_parameters==NULL)
	{
		int num_matches = read.length() ;
		
		for (size_t i=0; i<strlen(read_anno); i++)
		{
			if (read_anno[i]=='[')
			{
				num_matches-- ;
				i+=3 ;
			}
		}
		return num_matches ;
	}

	int num_qual_support = 0;
	if (alignment_parameters->num_qualityPlifs > 0)
		num_qual_support = alignment_parameters->qualityPlifs[0].len;
//	Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, true, verbosity);
	Alignment alignment(alignment_parameters->num_qualityPlifs, num_qual_support, !_config.NO_QPALMA, verbosity);

	double* prb = NULL ;
	try
	{
		prb = new double[read.length()] ;
	}
	catch (std::bad_alloc&)
	{
		fprintf(stderr,	"[score_unspliced] Could not allocate memory (_read.id() = %s)\n", read.id());
		return -1;
	} 
	int offset_switched=0;
	
	for (size_t i = 0; i < read.length(); i++)
	{
		prb[i] = (read.quality(0)[i] - alignment_parameters->quality_offset);
		if (prb[i]<-10 || prb[i]>70)
			fprintf(stderr, "prb[%i]=%f (offset=%i, %s, %s)\n", (int)i, prb[i], alignment_parameters->quality_offset, read.quality(0), read.data()) ;

		if (_config.QPALMA_PRB_OFFSET_FIX)
		{
			if (alignment_parameters->quality_offset==33 && prb[i]>70)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 64) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=64 ;
				read.set_quality_offset(64) ;
				offset_switched=1;
				
			}
			if (alignment_parameters->quality_offset==64 && prb[i]<-10)
			{
				fprintf(stderr, "setting prb offset from %i to %i\n", alignment_parameters->quality_offset, 33) ;
				_config.QPALMA_PRB_OFFSET_FIX=false ;
				alignment_parameters->quality_offset=33 ;
				read.set_quality_offset(33) ;
				offset_switched=1;
			}
		}

		//assert(prb[i]>=-10 && prb[i]<=70) ;
	}

	// If offset has been switched: recompute qualities for all positions
	if (offset_switched==1){
		for (size_t i = 0; i < read.length(); i++)
			prb[i] = (read.quality(0)[i] - alignment_parameters->quality_offset);
	}
	

	if (ori=='-')
		reverse(prb,read.length());

			
	double score1 = alignment.scoreUnsplicedAlignment(read_anno, prb, read.length(), alignment_parameters->qualityPlifs, alignment_parameters->matchmatrix,strand) ;
	//float score2 = alignment.scoreUnsplicedAlignment(read_anno, prb, read.length(), alignment_parameters->qualityPlifs, alignment_parameters->matchmatrix, '-') ;
	
	delete[] prb ;
	//if (score1>score2)
	return score1 ;
		//return score2 ;
}
