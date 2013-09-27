#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>

#include <palmapper/Config.h>
#include <palmapper/QueryFile.h>
#include <palmapper/Read.h>
#include <palmapper/Genome.h>

inline int max(int a, int b)
{
	if (a>b)
		return a ;
	else
		return b ;
}

// We need to access the name of the program before reaching main()
extern char const *__progname;
Personality Config::getPersonality() {
	return ::strcmp("genomemapper", __progname) == 0 ? GenomeMapper : Palmapper;
}


Config::Config() {
	_personality = getPersonality();

	BWA_INDEX=0; // Bwt	or genomemapper index
	NUM_THREADS = 1;//::sysconf(_SC_NPROCESSORS_ONLN);
	OUTPUT_FILTER = OUTPUT_FILTER_DEFAULT ;
	OUTPUT_FILTER_NUM_TOP = 10 ;
	OUTPUT_FILTER_NUM_LIMIT = 0 ; // all

	FIXTRIM_STRATEGY_LEN = 10000000;
	FIXTRIMRIGHT_STRATEGY_LEN = 0 ;
	FIXTRIMLEFT_STRATEGY_LEN = 0 ;
	RTRIM_STRATEGY=0 ;
	RTRIM_STRATEGY_MIN_LEN = 25 ;
	RTRIM_STRATEGY_STEP = DEFAULT_SETTING ;
	POLYTRIM_STRATEGY=0 ;
	POLYTRIM_STRATEGY_MIN_LEN=25 ;
	POLYTRIM_STRATEGY_STEP = DEFAULT_SETTING ;
	POLYTRIM_STRATEGY_POLY_MIN_LEN = 10 ;
	ADAPTERTRIM_STRATEGY = 0 ;
	ADAPTERTRIM_STRATEGY_MIN_LEN = 40 ;
	ADAPTERTRIM_STRATEGY_LOG = std::string("") ;

	//int SUMMARY_HIT_STRATEGY_NUM_EDIT_OPS[2] ;
	HITLEN_LIMIT = 0;
	VERBOSE = 0;
	MAP_REVERSE = 1 ;
	STRINGENT_GAPLIMIT = 1 ;
	PRINT_SEQ = 0;
	INDEX_DEPTH = 0;
	INDEX_DEPTH_EXTRA = 1 ;
	INDEX_DEPTH_EXTRA_THRESHOLD = 1000000000 ;
	SEED_HIT_CANCEL_THRESHOLD = 100000000 ;
	OUTPUT_FORMAT = OUTPUT_FORMAT_DEFAULT ;
	OUTPUT_FORMAT_FLAGS = OUTPUT_FORMAT_FLAGS_DEFAULT ;
	OUTPUT_FORMAT_OPTION = OUTPUT_FORMAT_OPTION_DEFAULT ;

	FIRST_READ_NR = 0;
	LAST_READ_NR = 10000000;

	REPORT_FILE = NULL;
	REPORT_FILE_READONLY = 0 ;
	REPORT_REPETITIVE_SEEDS = 0 ;
	REPORT_MAPPED_REGIONS = 1 ;
	REPORT_MAPPED_READS = 1 ;
	REPORT_SPLICED_READS = 1 ;
	REPORT_RESET = 0 ;
	REPORT_GENOME_COVERAGE = 0 ;
	REPORT_GENOME_COVERAGE_FILE = NULL;
	
	REPORT_JUNCTIONS_FILE = std::string("") ;
	REPORT_JUNCTIONS=0;
	MAP_JUNCTIONS_FILE = std::string("") ;
	MAP_JUNCTIONS=0;
	MAP_JUNCTIONS_COVERAGE=2;

	QPALMA_USE_MAP = 1 ;
	QPALMA_USE_MAP_MAX_SIZE = 10000 ;
	QPALMA_USE_SPLICE_SITES = 0 ;
	QPALMA_USE_SPLICE_SITES_THRESH_DON = 0.0 ;
	QPALMA_USE_SPLICE_SITES_THRESH_ACC = 0.0 ;
	QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC = 0.0 ;
	QPALMA_MIN_NUM_MATCHES = DEFAULT_SETTING ;
	QPALMA_PRB_OFFSET_FIX = false ;
	
	ANNOTATED_SPLICE_SITES_FILE= std::string("") ;
	SCORE_ANNOTATED_SPLICE_SITES = 0 ;

	READ_COUNT_LIMIT = 0 ; // limits the number of reads for alignment
	LOG_TRIGGERED = false;  // #A#
	FILTER_BY_MAX_MISMATCHES = 0 ;
	FILTER_BY_MAX_GAPS = 0 ;
	FILTER_BY_SPLICE_SITES = true ;
	FILTER_BY_SPLICE_SITES_REGION = 5 ;
	FILTER_BY_SPLICE_SITES_EDIT_MIN = 0 ;
	FILTER_BY_SPLICE_SITES_THRESH_ACC=0 ;
	FILTER_BY_SPLICE_SITES_THRESH_DON=0 ;
	FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC = 0.01;
	NO_SPLICE_PREDICTIONS=0 ;
	INDEX_PRECACHE = 0 ;
	FLANKING = 0;
	NUM_EDIT_OPS = DEFAULT_SETTING ;
	NUM_MISMATCHES = DEFAULT_SETTING ;
	NUM_GAPS = DEFAULT_SETTING ;
	MM_SCORE = 4;
	M_SCORE = 0;
	GAP_SCORE = 5;
	GAPS_MOST_RIGHT = 0;
	OVERHANG_ALIGNMENT = 0;
	SCORES_OUT = 1;
	SPLICED_HITS = 0 ;
	SPLICED_HIT_MIN_LENGTH_SHORT = DEFAULT_SETTING ;
	SPLICED_HIT_MIN_LENGTH_COMB = DEFAULT_SETTING ;
	SPLICED_HIT_MIN_LENGTH_LONG = DEFAULT_SETTING ;
	SPLICED_LONGEST_INTRON_LENGTH = DEFAULT_SETTING ;
	SPLICED_SHORTEST_INTRON_LENGTH = DEFAULT_SETTING ;
	SPLICED_MAX_NUM_ALIGNMENTS = 10 ;
	SPLICED_CLUSTER_TOLERANCE = 10 ;
	SPLICED_MAX_INTRONS = DEFAULT_SETTING ;
	SPLICED_MIN_SEGMENT_LENGTH = DEFAULT_SETTING ;

	STATISTICS = 0;
	CHROM_CONTAINER_SIZE = 15000000 ;

	BEST_HIT_STRATEGY = 0 ;
	ALL_HIT_STRATEGY = 0 ;
	SUMMARY_HIT_STRATEGY= 0 ;
	HITLEN_LIMIT = 0 ;

	BSSEQ = 0;

	LEFTOVER_FILE_NAME = std::string("/dev/null") ;
	SAMTOOLS_PATH_NAME = std::string("") ;

	ACC_CONSENSUS.push_back(strdup("AG")) ;
	ACC_CONSENSUS_REV.push_back(strdup("CT")) ;
	DON_CONSENSUS.push_back(strdup("GT")) ;
	DON_CONSENSUS.push_back(strdup("GC")) ;
	DON_CONSENSUS_REV.push_back(strdup("AC")) ;
	DON_CONSENSUS_REV.push_back(strdup("GC")) ;

	non_consensus_search = false ;
	non_consensus_search_gap=1 ;
	non_consensus_search_discount=1 ;
	// Number of additional matches you have to find in case of a non consensus search compared to consensus one (based on QMM)
	MIN_NUM_MATCHES_PEN=2;
	
	INCLUDE_UNMAPPED_READS_SAM=false;
	
	NO_QPALMA = false;
	
    STRAND = -1 ;
	PROTOCOL = PROTOCOL_UNSTRANDED ;

	Q_QUERY_FILE_NAMES= std::string("") ;
	Q1_QUERY_FILE_NAMES= std::string("") ;
	Q2_QUERY_FILE_NAMES= std::string("") ;

};

int Config::applyDefaults(Genome * genome)
{
	int any_default=false ;
	
	if ((SPLICED_HITS && (SPLICED_HIT_MIN_LENGTH_SHORT == DEFAULT_SETTING || SPLICED_HIT_MIN_LENGTH_LONG == DEFAULT_SETTING || 
						  SPLICED_HIT_MIN_LENGTH_COMB == DEFAULT_SETTING || SPLICED_MAX_INTRONS == DEFAULT_SETTING)) || 
		NUM_EDIT_OPS == DEFAULT_SETTING || NUM_MISMATCHES == DEFAULT_SETTING || NUM_GAPS == DEFAULT_SETTING || OUTPUT_FILTER == OUTPUT_FILTER_DEFAULT ||
		(SPLICED_HITS && (SPLICED_LONGEST_INTRON_LENGTH == DEFAULT_SETTING || SPLICED_SHORTEST_INTRON_LENGTH == DEFAULT_SETTING || 
						  SPLICED_MIN_SEGMENT_LENGTH==DEFAULT_SETTING)) || 
		(OUTPUT_FORMAT==OUTPUT_FORMAT_DEFAULT) ||
		((int)POLYTRIM_STRATEGY_STEP == DEFAULT_SETTING && POLYTRIM_STRATEGY) || 
		((int)RTRIM_STRATEGY_STEP == DEFAULT_SETTING && RTRIM_STRATEGY))
	{
		fprintf(stdout, "\nSetting default parameters:\n") ;
		any_default=true ;
	}

	if (_personality == Palmapper)  {
		int read_length = QueryFile::determine_read_length(QUERY_FILE_NAMES,QUERY_FILE_STRANDS);
		{
			bool line_started=false ;
			if ((SPLICED_HITS && (SPLICED_HIT_MIN_LENGTH_SHORT == DEFAULT_SETTING || SPLICED_HIT_MIN_LENGTH_LONG == DEFAULT_SETTING || SPLICED_HIT_MIN_LENGTH_COMB == DEFAULT_SETTING || SPLICED_MAX_INTRONS == DEFAULT_SETTING)) ||
				NUM_EDIT_OPS == DEFAULT_SETTING || NUM_MISMATCHES == DEFAULT_SETTING || NUM_GAPS == DEFAULT_SETTING)
			{
				line_started=true ;
				fprintf(stdout, "* Automatically determining alignment parameters based on read length (%int):", read_length) ;
			}
			if (SPLICED_HITS && SPLICED_HIT_MIN_LENGTH_SHORT == DEFAULT_SETTING)
			{
				SPLICED_HIT_MIN_LENGTH_SHORT = 15 ;
				fprintf(stdout, " -K %i", SPLICED_HIT_MIN_LENGTH_SHORT) ;
			}
			if (SPLICED_HITS && SPLICED_HIT_MIN_LENGTH_LONG == DEFAULT_SETTING)
			{
				SPLICED_HIT_MIN_LENGTH_LONG = (read_length/4<20) ? 20 : read_length/4 ;
				fprintf(stdout, " -L %i", SPLICED_HIT_MIN_LENGTH_LONG) ;
			}
			if (SPLICED_HITS && SPLICED_HIT_MIN_LENGTH_COMB == DEFAULT_SETTING)
			{
				SPLICED_HIT_MIN_LENGTH_COMB = (read_length/2<30) ? 30 : read_length/2 ;
				if (SPLICED_HIT_MIN_LENGTH_COMB<SPLICED_HIT_MIN_LENGTH_LONG)
					SPLICED_HIT_MIN_LENGTH_COMB=SPLICED_HIT_MIN_LENGTH_LONG ;
				fprintf(stdout, " -C %i", SPLICED_HIT_MIN_LENGTH_COMB) ;
			}
			if (SPLICED_HITS && SPLICED_MAX_INTRONS == DEFAULT_SETTING)
			{
				SPLICED_MAX_INTRONS = (read_length>50) ? (read_length>=100 ? 3 : 2) : 1 ;
				if (read_length>200) SPLICED_MAX_INTRONS = 6 ;
				fprintf(stdout, " -NI %i", SPLICED_MAX_INTRONS) ;
			}
			if (NUM_EDIT_OPS == DEFAULT_SETTING)
			{
				NUM_EDIT_OPS = read_length*0.08 ;
				fprintf(stdout, " -E %i", NUM_EDIT_OPS) ;
			}
			if (NUM_MISMATCHES == DEFAULT_SETTING)
			{
				NUM_MISMATCHES = read_length*0.08 ;
				fprintf(stdout, " -M %i", NUM_MISMATCHES) ;
			}
			if (NUM_GAPS == DEFAULT_SETTING)
			{
				NUM_GAPS = read_length*0.03 ;
				fprintf(stdout, " -G %i", NUM_GAPS) ;
			}
			if (line_started)
				fprintf(stdout, "\n") ;

			if (SPLICED_HITS && QPALMA_MIN_NUM_MATCHES == DEFAULT_SETTING)
			{
				QPALMA_MIN_NUM_MATCHES=5 ;
				//if (non_consensus_search)
				//	QPALMA_MIN_NUM_MATCHES+=2 ;
				fprintf(stdout, "* Automatically determined minimal match length near splice sites (%int)\n", QPALMA_MIN_NUM_MATCHES) ;
			}

			if (SPLICED_HITS && SPLICED_MIN_SEGMENT_LENGTH == DEFAULT_SETTING)
			{
				SPLICED_MIN_SEGMENT_LENGTH=QPALMA_MIN_NUM_MATCHES ;

				if (read_length>40)
					SPLICED_MIN_SEGMENT_LENGTH = max(QPALMA_MIN_NUM_MATCHES, 6) ;
				if (read_length>=75)
					SPLICED_MIN_SEGMENT_LENGTH = max(QPALMA_MIN_NUM_MATCHES, 8) ;
				if (read_length>=100)
					SPLICED_MIN_SEGMENT_LENGTH = max(QPALMA_MIN_NUM_MATCHES, 10) ;

				fprintf(stdout, "* Automatically determined minimal segment length in spliced alignments based on read length (%int)\n",
						SPLICED_MIN_SEGMENT_LENGTH) ;
			}
		}

		if (OUTPUT_FILTER == OUTPUT_FILTER_DEFAULT)
		{
			OUTPUT_FILTER = OUTPUT_FILTER_TOP ;
			OUTPUT_FILTER_NUM_TOP = 5 ;
			fprintf(stdout, "* Reporting the best %i alignments per read\n", OUTPUT_FILTER_NUM_TOP) ;
		}

		if (SPLICED_HITS && SPLICED_LONGEST_INTRON_LENGTH == DEFAULT_SETTING)
		{
			unsigned long int genome_size = 0 ;
			for (unsigned int i=0; i<genome->nrChromosomes(); i++)
				genome_size+=genome->chromosome(i).length() ;

			if (genome_size<900000000)
				SPLICED_LONGEST_INTRON_LENGTH = 20000 ;
			else if (genome_size<500000000)
				SPLICED_LONGEST_INTRON_LENGTH = 50000 ;
			else
				SPLICED_LONGEST_INTRON_LENGTH = 200000 ;
			fprintf(stdout, "* Automatically determined maximal intron size based on genome size (%ikb)\n", SPLICED_LONGEST_INTRON_LENGTH/1000) ;
		}

		if (SPLICED_HITS && SPLICED_SHORTEST_INTRON_LENGTH == DEFAULT_SETTING)
		{
			SPLICED_SHORTEST_INTRON_LENGTH = 30 ;
			fprintf(stdout, "* Automatically determined minimal intron size(%int)\n", SPLICED_SHORTEST_INTRON_LENGTH) ;
		}

		// determine default output format
		if (OUTPUT_FORMAT==OUTPUT_FORMAT_DEFAULT) {
			if (SPLICED_HITS)
			{
				OUTPUT_FORMAT=OUTPUT_FORMAT_SAM ;
				fprintf(stdout, "* Selecting SAM output format\n") ;
			}
			else
			{
				OUTPUT_FORMAT=OUTPUT_FORMAT_SAM ;
				fprintf(stdout, "* Selecting SAM output format\n") ;
			}
		}
	} else {

		if (OUTPUT_FORMAT==OUTPUT_FORMAT_DEFAULT) OUTPUT_FORMAT = OUTPUT_FORMAT_SHORE;

		if (NUM_EDIT_OPS == DEFAULT_SETTING) NUM_EDIT_OPS = 3;
		if (NUM_MISMATCHES == DEFAULT_SETTING) NUM_MISMATCHES = 3;
		if (NUM_GAPS == DEFAULT_SETTING) NUM_GAPS = 1;

		if (OUTPUT_FILTER == OUTPUT_FILTER_DEFAULT) {
			BEST_HIT_STRATEGY = 1;
		}
	}

	if (POLYTRIM_STRATEGY)
	{
		if ((int)POLYTRIM_STRATEGY_STEP == DEFAULT_SETTING)
		{
			POLYTRIM_STRATEGY_STEP = (NUM_EDIT_OPS>=1) ? NUM_EDIT_OPS : 1 ;
			fprintf(stdout, "* Automatically selecting polytrim step size: %int\n", POLYTRIM_STRATEGY_STEP) ;
		}
	}

	if (RTRIM_STRATEGY)
	{
		if ((int)RTRIM_STRATEGY_STEP == DEFAULT_SETTING)
		{
			RTRIM_STRATEGY_STEP = (NUM_EDIT_OPS>=2) ? NUM_EDIT_OPS/2 : 1 ;
			fprintf(stdout, "* Automatically selecting rtrim step size: %int\n", RTRIM_STRATEGY_STEP) ;
		}
	}

	if (any_default)
		fprintf(stdout, "\n") ;
	

	return 0 ;
}

int Config::checkConfig()
{
	if (_personality == Palmapper) {
		int read_length = QueryFile::determine_read_length(QUERY_FILE_NAMES,QUERY_FILE_STRANDS);

		if (SPLICED_OUT_FILE_NAME.length()>0 && !SPLICED_HITS)
		{
			fprintf(stderr, "ERROR: output files for spliced hits provided, but no spliced alignment is performed\n");
			exit(1);
		}

		if (SPLICED_HITS && !(NO_SPLICE_PREDICTIONS || (ACC_FILES.length()>0 && DON_FILES.length()>0) || ANNOTATED_SPLICE_SITES_FILE.length()>0))
		{
			fprintf(stderr, "ERROR: for spliced alignments either -acc and -don and/or -score-annotated-splice-sites or -no-ss-pred need to be given as argument\n");
			exit(1);
		}

		if (NO_SPLICE_PREDICTIONS && (ACC_FILES.length()>0 || DON_FILES.length()>0 || ANNOTATED_SPLICE_SITES_FILE.length()>0))
		{
			fprintf(stderr, "ERROR: the options -acc/-don/-score-annotated-splice-sites and -no-ss-pred have to be used exclusively\n");
			exit(1);
		}

		if (SPLICED_HITS && (OUTPUT_FORMAT==OUTPUT_FORMAT_SHORE || OUTPUT_FORMAT==OUTPUT_FORMAT_BED))
		{
			fprintf(stderr, "ERROR: SHORE or BED format currently do not support spliced alignments (choose BEDX or SAM) %i\n", (int)OUTPUT_FORMAT) ;
			exit(1) ;
		}

		if(SPLICED_HITS && (SPLICED_HIT_MIN_LENGTH_LONG>read_length/2)){
			fprintf(stderr,"WARNING: Minimal length of long hit (%i) is greater to the half of read length. Reset to half read length (%i)\n", SPLICED_HIT_MIN_LENGTH_LONG, read_length);
			SPLICED_HIT_MIN_LENGTH_LONG=read_length/2;
		}
	}

	if (RTRIM_STRATEGY && POLYTRIM_STRATEGY)
	{
		fprintf(stderr, "ERROR: RTRIM and POLYTRIM cannot be combined\n") ;
		exit(1) ;
	}


	if (INCLUDE_UNMAPPED_READS_SAM && OUTPUT_FORMAT!=OUTPUT_FORMAT_SAM && OUTPUT_FORMAT!=OUTPUT_FORMAT_BAM)
	{
		fprintf(stderr, "WARNING: unmapped reads can only be written in the same output file than mapped reads for sam format\n") ; 
		INCLUDE_UNMAPPED_READS_SAM=false;
		
	}

	if (MAX_EDIT_OPS<NUM_EDIT_OPS)
	{
		fprintf(stderr, "ERROR: NUM_EDIT_OPS>MAX_EDIT_OPS, please increase MAX_EDIT_OPS in config.h\n") ;
		exit(1) ;
	}
	
	if (SPLICED_MAX_INTRONS>1 && non_consensus_search)
	{
		fprintf(stdout, "WARNING: non-consensus search is slow with more than one intron per alignment\n") ;
	}

	if (FIXTRIM_STRATEGY_LEN< 10000000 && (FIXTRIMLEFT_STRATEGY_LEN>0 || FIXTRIMRIGHT_STRATEGY_LEN>0))
	{
		fprintf(stderr,	"ERROR: -fixtrim and -fixtrim[right|left] options cannot be combined\n");
		exit(1) ;
	}
	
	return 0 ;
}

int Config::postprocess_query_filenames(std::string filenames, int strand) 
{

	if (filenames.length()==0)
		return 0;
		
	int previousfound=0;
	int found=filenames.find(",");
	std::string filename;
	
	while (found >= 0)
	{
		
		QUERY_FILE_NAMES.push_back(filenames.substr(previousfound, found-previousfound));
		QUERY_FILE_STRANDS.push_back(strand);
		
		previousfound=found+1;
		found=filenames.find(",",found+1);
	}
		
	QUERY_FILE_NAMES.push_back(filenames.substr(previousfound));
	QUERY_FILE_STRANDS.push_back(strand);
	
	return 0;
	
}

int Config::postprocess_consensus_list(std::vector<const char *> & consensus_list) 
{
  for (unsigned int i=0; i<consensus_list.size(); i++)
    {
      if (strlen(consensus_list[i])!=2)
	{
	  fprintf(stderr, "invalid consensus %s\n", consensus_list[i]) ;
	  return -1 ;
	}
      if (consensus_list[i][0]=='?')
	{
	  char con[3] ;
	  con[2]=0 ;
	  con[1]=consensus_list[i][1] ;
	  const char *acgt="ACGT" ;
	  con[0]=acgt[0] ;
	  //free(consensus_list[i]) ;
	  consensus_list[i]=strdup(con) ;
	  for (int j=1; j<4; j++)
	    {
	      con[0]=acgt[j] ;
	      consensus_list.push_back(strdup(con)) ;
	    }
	}
      if (consensus_list[i][1]=='?')
	{
	  char con[3] ;
	  con[2]=0 ;
	  con[0]=consensus_list[i][0] ;
	  const char *acgt="ACGT" ;
	  con[1]=acgt[0] ;
	  //free(consensus_list[i]) ;
	  consensus_list[i]=strdup(con) ;
	  for (int j=1; j<4; j++)
	    {
	      con[1]=acgt[j] ;
	      consensus_list.push_back(strdup(con)) ;
	    }
	}
    }

  fprintf(stdout, "consensus list: ") ;
  for (unsigned int i=0; i<consensus_list.size(); i++)
    fprintf(stdout, "%s, ", consensus_list[i]) ;

  return 0 ;
}

int Config::parseCommandLine(int argc, char *argv[]) 
{
	int i;
	char not_defined;
	char has_index = 0;
	char has_genome = 0;

	for (i = 1; i < argc; i++) {
		not_defined = 1;

		//genome file
		if (strcmp(argv[i], "-i") == 0) {

			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -i\n") ;
				usage();
				exit(1);
			}
			i++;
			GENOME_FILE_NAME.assign(argv[i]);
			has_genome = 1;
			has_index = 1;

			CHR_INDEX_FILE_NAME.assign(argv[i]) ;
			CHR_INDEX_FILE_NAME += ".cid";

			META_INDEX_FILE_NAME.assign(argv[i]);
			META_INDEX_FILE_NAME += ".mta";

			INDEX_FWD_FILE_NAME.assign(argv[i]);
			INDEX_FWD_FILE_NAME += ".mfd";

			//INDEX_REV_FILE_NAME.assign(argv[i]);
			//INDEX_REV_FILE_NAME += ".mrc";

		}

		/**

		 //chr index file
		 if(strcmp(argv[i],"-x")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) 
		 { 
		 fprintf(stderr, "ERROR: Argument missing for option -i\n") ;
		 usage(); exit(1); 
		 }
		 i++;
		 strcpy(CHR_INDEX_FILE_NAME, argv[i]);
		 has_index = 1;
		 }

		 //meta index file
		 if(strcmp(argv[i],"-t")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) { usage(); exit(1); }
		 i++;
		 strcpy(META_INDEX_FILE_NAME, argv[i]);
		 has_meta_index = 1;
		 }

		 //index fwd file
		 if(strcmp(argv[i],"-z")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) { usage(); exit(1); }
		 i++;
		 strcpy(INDEX_FWD_FILE_NAME, argv[i]);
		 has_fwd_index = 1;
		 }

		 //index rev file
		 if(strcmp(argv[i],"-y")==0){
		 not_defined = 0;
		 if(i+1 > argc - 1) { usage(); exit(1); }
		 i++;
		 strcpy(INDEX_REV_FILE_NAME, argv[i]);
		 has_rev_index = 1;
		 }

		 */

		//BWA index
		if (strcmp(argv[i], "-bwa") == 0) {
			not_defined = 0;
			BWA_INDEX = 1 ;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -bwa\n") ;
				usage();
				exit(1);
			}
			i++;
			INDEX_DEPTH = atoi(argv[i]) ;

		}
		

		 //nr threads
		char const *threadOpt = _personality == Palmapper ? "-threads" : "-t";
		if(strcmp(argv[i],threadOpt)==0) {
			not_defined = 0;
			if(i+1 > argc - 1) 
			{
				fprintf(stderr, "ERROR: Argument missing for option %s\n", threadOpt) ;
				usage(); 
				exit(1); 
			}
			i++;
		    NUM_THREADS = atoi(argv[i]) ;
			if (NUM_THREADS<1)
			{
				fprintf(stderr,	"ERROR: number of threads too small\n");
				exit(1) ;
			}
		}

		//query file
		if (strcmp(argv[i], "-q") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -q\n") ;
				usage();
				exit(1);
			}
			i++;
			Q_QUERY_FILE_NAMES=strdup(argv[i]);
		}

		//output file
		if (strcmp(argv[i], "-o") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -o\n") ;
				usage();
				exit(1);
			}
			i++;
			OUT_FILE_NAME.assign(argv[i]);
		}

		//fix read trimming
		if (strcmp(argv[i], "-fixtrim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -fixtrim\n") ;
				usage();
				exit(1);
			}
			i++;
			FIXTRIM_STRATEGY_LEN = atoi(argv[i]) ;
			if (FIXTRIM_STRATEGY_LEN<INDEX_DEPTH)
			{
				fprintf(stderr,	"ERROR: minimal fixtrim alignment length too short\n");
				exit(1) ;
			}
		}


		//fix read trimming left
		if (strcmp(argv[i], "-fixtrimleft") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -fixtrimleft\n") ;
				usage();
				exit(1);
			}
			i++;
			FIXTRIMLEFT_STRATEGY_LEN = atoi(argv[i]) ;
		}

		//fix read trimming right
		if (strcmp(argv[i], "-fixtrimright") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -fixtrimright\n") ;
				usage();
				exit(1);
			}
			i++;
			FIXTRIMRIGHT_STRATEGY_LEN = atoi(argv[i]) ;
		}

		//partial alignments for rtrim
		if (strcmp(argv[i], "-rtrim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -rtrim\n") ;
				usage();
				exit(1);
			}
			i++;
			RTRIM_STRATEGY = 1 ;
			RTRIM_STRATEGY_MIN_LEN = atoi(argv[i]) ;
			if (RTRIM_STRATEGY_MIN_LEN<INDEX_DEPTH)
			{
				fprintf(stderr,	"ERROR: minimal rtrim alignment length too short\n");
				exit(1) ;
			}
		}

		//partial alignments for rtrim
		if (strcmp(argv[i], "-rtrim-step") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -rtrim-step\n") ;
				usage();
				exit(1);
			}
			i++;
			RTRIM_STRATEGY = 1 ;
			RTRIM_STRATEGY_STEP = atoi(argv[i]) ;
			if (RTRIM_STRATEGY_STEP<1)
			{
				fprintf(stderr,	"ERROR: rtrim step size too short\n");
				exit(1) ;
			}
		}

		//polyA/T trimming
		if (strcmp(argv[i], "-polytrim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -polytrim\n") ;
				usage();
				exit(1);
			}
			i++;
			POLYTRIM_STRATEGY = 1 ;
			POLYTRIM_STRATEGY_MIN_LEN = atoi(argv[i]) ;
			if (POLYTRIM_STRATEGY_MIN_LEN<INDEX_DEPTH)
			{
				fprintf(stderr,	"ERROR: minimal polytrim alignment length too short\n");
				exit(1) ;
			}
		}

		//adapter trimming
		if (strcmp(argv[i], "-adaptertrim") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -adaptertrim\n") ;
				usage();
				exit(1);
			}
			i++;
			ADAPTERTRIM_STRATEGY = 1 ;
			ADAPTERTRIM_STRATEGY_MIN_LEN = atoi(argv[i]) ;
			if (ADAPTERTRIM_STRATEGY_MIN_LEN<INDEX_DEPTH)
			{
				fprintf(stderr,	"ERROR: minimal polytrim alignment length too short\n");
				exit(1) ;
			}
		}

		//adapter trimming
		if (strcmp(argv[i], "-adaptertrim-log") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -adaptertrim-log\n") ;
				usage();
				exit(1);
			}
			i++;
			ADAPTERTRIM_STRATEGY_LOG = std::string(argv[i]) ;
		}

		if (_personality == Palmapper) {
			//report output file
			if (strcmp(argv[i], "-report") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report\n") ;
					usage();
					exit(1);
				}
				i++;
				REPORT_FILE=strdup(argv[i]) ;
				REPORT_FILE_READONLY = 0 ;
			}

			//report output file
			if (strcmp(argv[i], "-report-ro") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report-ro\n") ;
					usage();
					exit(1);
				}
				i++;
				REPORT_FILE=strdup(argv[i]) ;
				REPORT_FILE_READONLY = 1 ;
			}

			//report output file
			if (strcmp(argv[i], "-report-coverage-wig") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report-coverage-wig\n") ;
					usage();
					exit(1);
				}
				i++;
				REPORT_GENOME_COVERAGE_FILE=strdup(argv[i]) ;
				REPORT_MAPPED_READS = 1 ;
				REPORT_SPLICED_READS = 1 ;
				REPORT_GENOME_COVERAGE = 2 ;
			}
			//report output file
			if (strcmp(argv[i], "-report-coverage-map") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report-coverage-map\n") ;
					usage();
					exit(1);
				}
				i++;
				REPORT_GENOME_COVERAGE_FILE=strdup(argv[i]) ;
				REPORT_MAPPED_READS = 1 ;
				REPORT_SPLICED_READS = 1 ;
				REPORT_GENOME_COVERAGE = 1 ;
			}

			if (strcmp(argv[i], "-score-annotated-splice-sites") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -score-annotated-splice-sites\n") ;
					usage();
					exit(1);
				}
				i++;
				ANNOTATED_SPLICE_SITES_FILE=strdup(argv[i]) ;
				SCORE_ANNOTATED_SPLICE_SITES = 1 ;
			}

			//report junctions
			if (strcmp(argv[i], "-report-junctions") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report-junctions\n") ;
					usage();
					exit(1);
				}
				i++;
				REPORT_JUNCTIONS_FILE=strdup(argv[i]) ;
				REPORT_JUNCTIONS = 1 ;
			}
			if (strcmp(argv[i], "-junction-remapping") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -junction-remapping\n") ;
					usage();
					exit(1);
				}
				i++;
				MAP_JUNCTIONS_FILE=strdup(argv[i]) ;
				MAP_JUNCTIONS = 1 ;
			}


			if (strcmp(argv[i], "-junction-remapping-coverage") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -junction-remapping-coverage\n") ;
					usage();
					exit(1);
				}
				i++;
				MAP_JUNCTIONS_COVERAGE=atoi(argv[i]) ;
				assert(MAP_JUNCTIONS_COVERAGE>=0) ;
			}

			//report repetitive seeds
			if (strcmp(argv[i], "-report-rep-seed") == 0) {
				not_defined = 0;
				REPORT_REPETITIVE_SEEDS = 1 ;
				//assert(REPORT_SPLICE_SITES==0) ; // currently not supported
			}
			if (strcmp(argv[i], "-no-report-rep-seed") == 0) {
				not_defined = 0;
				REPORT_REPETITIVE_SEEDS = 0 ;
				//assert(REPORT_SPLICE_SITES==0) ; // currently not supported
			}

			//report mapped regions
			if (strcmp(argv[i], "-report-map-region") == 0) {
				not_defined = 0;
				REPORT_MAPPED_REGIONS  = 1 ;
			}
			if (strcmp(argv[i], "-no-report-map-region") == 0) {
				not_defined = 0;
				REPORT_MAPPED_REGIONS  = 0 ;
			}

			//report mapped regions
			if (strcmp(argv[i], "-report-map-read") == 0) {
				not_defined = 0;
				REPORT_MAPPED_READS = 1 ;
			}
			if (strcmp(argv[i], "-no-report-map-read") == 0) {
				not_defined = 0;
				REPORT_MAPPED_READS = 0 ;
			}

			//report mapped regions
			if (strcmp(argv[i], "-report-spliced-read") == 0) {
				not_defined = 0;
				REPORT_SPLICED_READS = 1 ;
			}
			if (strcmp(argv[i], "-no-report-spliced-read") == 0) {
				not_defined = 0;
				REPORT_SPLICED_READS = 0 ;
			}

			//report splice sites - confidence threshold
			if (strcmp(argv[i], "-report-splice-sites") == 0) {
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report-splice-sites\n") ;
					usage();
					exit(1);
				}
				i++;
				QPALMA_USE_SPLICE_SITES_THRESH_ACC = atof(argv[i]);
				QPALMA_USE_SPLICE_SITES_THRESH_DON = atof(argv[i]);
				QPALMA_USE_SPLICE_SITES= 1 ;
				not_defined = 0;
				//assert(REPORT_REPETITIVE_SEEDS==0) ; // currently not supported
				//assert(REPORT_SPLICE_SITES_THRESH_TOP_PERC==0.0) ;
			}

			//report splice sites - percentile threshold
			if (strcmp(argv[i], "-report-splice-sites-top-perc") == 0) {
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report-splice-sites-top-perc\n") ;
					usage();
					exit(1);
				}
				i++;
				QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC = atof(argv[i]);
				assert(QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC>=0 && QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC<=1.0) ;
				QPALMA_USE_SPLICE_SITES= 1 ;
				not_defined = 0;
				//assert(REPORT_REPETITIVE_SEEDS==0) ; // currently not supported
				//assert(REPORT_SPLICE_SITES_THRESH==0.0) ;
			}

			// filter by splice sites - percentile threshold
			if (strcmp(argv[i], "-filter-splice-sites-top-perc") == 0) {
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -filter-splice-sites-top-perc\n") ;
					usage();
					exit(1);
				}
				i++;
				FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC = atof(argv[i]);
				assert(FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC>=0 && FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC<=1.0) ;
				not_defined = 0;
				FILTER_BY_SPLICE_SITES_THRESH_ACC = 0.0 ;
				FILTER_BY_SPLICE_SITES_THRESH_DON = 0.0 ;
			}

			// filter by mismatches
			if (strcmp(argv[i], "-filter-max-mismatches") == 0) {
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -filter-max-mismatches\n") ;
					usage();
					exit(1);
				}
				i++;
				FILTER_BY_MAX_MISMATCHES = atoi(argv[i]);
				assert(FILTER_BY_MAX_MISMATCHES>=0) ;
				not_defined = 0;
			}

			// splice-site based filter: require at least this many edits
			if (strcmp(argv[i], "-filter-max-edit") == 0) {
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -filter-splice-min-edit\n") ;
					usage();
					exit(1);
				}
				i++;
				FILTER_BY_SPLICE_SITES_EDIT_MIN = atoi(argv[i]);
				assert(FILTER_BY_SPLICE_SITES_EDIT_MIN>=0) ;
				not_defined = 0;
			}

			// splice-site based filter: consider a region of this length around the read; -1 switches off the splice site-based filter
			if (strcmp(argv[i], "-filter-splice-region") == 0)
			{
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -filter-splice-region\n") ;
					usage();
					exit(1);
				}
				i++;
				FILTER_BY_SPLICE_SITES_REGION = atoi(argv[i]);
				if (FILTER_BY_SPLICE_SITES_REGION==-1)
					FILTER_BY_SPLICE_SITES = 0 ;

				assert(FILTER_BY_SPLICE_SITES_REGION>=-1) ;
				not_defined = 0;
			}

			// filter by mismatches
			if (strcmp(argv[i], "-filter-max-gaps") == 0) {
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -filter-max-gaps\n") ;
					usage();
					exit(1);
				}
				i++;
				FILTER_BY_MAX_GAPS = atoi(argv[i]);
				assert(FILTER_BY_MAX_GAPS>=0) ;
				not_defined = 0;
			}

			// reset the report (i.e. don't load it from file, even when available)
			if (strcmp(argv[i], "-report-reset") == 0) {
				not_defined = 0;
				REPORT_RESET = 1 ;
			}

			// use regions around mapped reads for qpalma alignment
			/*if (strcmp(argv[i], "-qpalma-use-map") == 0) {
				not_defined = 0;
				QPALMA_USE_MAP = 1 ;
				}*/

			// use regions around mapped reads for qpalma alignment
			if (strcmp(argv[i], "-qpalma-use-map-max-len") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -qpalma-use-map-max-len\n") ;
					usage();
					exit(1);
				}
				i++;
				QPALMA_USE_MAP_MAX_SIZE = atoi(argv[i]);
			}

			//spliced output file
			if (strcmp(argv[i], "-report-gff-init") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -report-gff-init\n") ;
					usage();
					exit(1);
				}
				i++;
				REPORT_GFF_FILE_NAME.assign(argv[i]);
			}

			//spliced output file
			if (strcmp(argv[i], "-H") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -H\n") ;
					usage();
					exit(1);
				}
				i++;
				SPLICED_OUT_FILE_NAME.assign(argv[i]);
			}

			if (strcmp(argv[i], "-read-id-prefix") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -read-id-prefix\n") ;
					usage();
					exit(1);
				}
				i++;
				READ_ID_PREFIX.assign(argv[i]);
			}

			//spliced hits min combined length
			if (strcmp(argv[i], "-C") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -C\n") ;
					usage();
					exit(1);
				}
				i++;
				SPLICED_HIT_MIN_LENGTH_COMB = atoi(argv[i]);
			}

			//spliced hits min length
			if (strcmp(argv[i], "-K") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -K\n") ;
					usage();
					exit(1);
				}
				i++;
				SPLICED_HIT_MIN_LENGTH_SHORT = atoi(argv[i]);
			}

			//spliced hits min length of longer hit
			if (strcmp(argv[i], "-L") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -L\n") ;
					usage();
					exit(1);
				}
				i++;
				SPLICED_HIT_MIN_LENGTH_LONG = atoi(argv[i]);
			}

			// longest intron length
			if (strcmp(argv[i], "-I") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -I\n") ;
					usage();
					exit(1);
				}
				i++;
				SPLICED_LONGEST_INTRON_LENGTH = atoi(argv[i]);
			}

			// longest intron length
			if (strcmp(argv[i], "-MI") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -MI\n") ;
					usage();
					exit(1);
				}
				i++;
				SPLICED_SHORTEST_INTRON_LENGTH = atoi(argv[i]);
			}

			// maximal number of spliced alignments to be performed per read
			if (strcmp(argv[i], "-SA") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -SA\n") ;
					usage();
					exit(1);
				}
				i++;
				int tmp = atoi(argv[i]);
				if (tmp < 0) {
					fprintf(stderr, "ERROR: Argument for option -SA too small\n") ;
					usage();
					exit(1);
				}
				SPLICED_MAX_NUM_ALIGNMENTS = tmp;
			}

			// maximal number of introns in spliced alignments
			if (strcmp(argv[i], "-NI") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -NI\n") ;
					usage();
					exit(1);
				}
				i++;
				int tmp = atoi(argv[i]);
				if (tmp < 0) {
					fprintf(stderr, "ERROR: Argument for option -NI too small\n") ;
					usage();
					exit(1);
				}
				SPLICED_MAX_INTRONS = tmp;
			}

			// How many matches are necessary for identifying a possible splice site in QPALMA recursive alignment algorithm?
			if (strcmp(argv[i], "-QMM") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -QMM\n") ;
					usage();
					exit(1);
				}
				i++;
				int tmp = atoi(argv[i]);
				if (tmp < 0) {
					fprintf(stderr, "ERROR: Argument for option -QMM too small\n") ;
					usage();
					exit(1);
				}
				QPALMA_MIN_NUM_MATCHES = tmp;
			}

			// how much distance to tolerate between a hit and an existing
			// hit cluster
			if (strcmp(argv[i], "-CT") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -CT\n") ;
					usage();
					exit(1);
				}
				i++;
				int tmp = atoi(argv[i]);
				if (tmp < 0) {
					fprintf(stderr, "ERROR: Argument for option -CT too small\n") ;
					usage();
					exit(1);
				}
				SPLICED_CLUSTER_TOLERANCE = tmp;
			}

			// limit the number of reads for alignment
			if (strcmp(argv[i], "-rlim") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -rlim\n") ;
					usage();
					exit(1);
				}
				i++;
				int tmp = atoi(argv[i]);
				if (tmp < 0) {
					fprintf(stderr, "ERROR: Argument for option -rlim too small\n") ;
					usage();
					exit(1);
				}
				READ_COUNT_LIMIT = tmp;
			}

			// limit the number of reads for alignment
			if (strcmp(argv[i], "-index-precache") == 0) {
				not_defined = 0;
				INDEX_PRECACHE = 1;
			}
		}

		//genome file
		if (strcmp(argv[i], "-samtools") == 0) {

			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -i\n") ;
				usage();
				exit(1);
			}
			i++;
			SAMTOOLS_PATH_NAME.assign(argv[i]);
			SAMTOOLS_PATH_NAME+="/";
			
		}

		//output format
		if (strcmp(argv[i], "-f") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -f\n") ;
				usage();
				exit(1);
			}
			std::string output = argv[++i];
			std::transform(output.begin(), output.end(), output.begin(), toupper);
			 
			if (output=="SHORE") {
				OUTPUT_FORMAT = OUTPUT_FORMAT_SHORE ;
			}
			else if (_personality == Palmapper && (output=="BED")) {
				OUTPUT_FORMAT = OUTPUT_FORMAT_BED;
			}
			else if (_personality == Palmapper && (output=="BEDX")) {
				OUTPUT_FORMAT = OUTPUT_FORMAT_BEDX;
			}
			else if (output=="SAM") {
				OUTPUT_FORMAT = OUTPUT_FORMAT_SAM;
			}
			else if (output=="BAM") {
				OUTPUT_FORMAT = OUTPUT_FORMAT_BAM;
			}
			else if (output=="BAMN") {
				OUTPUT_FORMAT = OUTPUT_FORMAT_BAM;
				OUTPUT_FORMAT_OPTION = OUTPUT_FORMAT_OPTION_SORTNAME ;
			}
			else if (output=="BAMP") {
				OUTPUT_FORMAT = OUTPUT_FORMAT_BAM;
				OUTPUT_FORMAT_OPTION = OUTPUT_FORMAT_OPTION_SORTPOS ;
			}
			else {
				fprintf(stderr,	"ERROR: Output file format %s not supported\n", output.c_str());
				exit(1);
			}
		}

		if (strcmp(argv[i], "-include-unmapped-reads") == 0) {
			not_defined = 0;
			INCLUDE_UNMAPPED_READS_SAM = true ;
		}

		//output format flags
		if (strcmp(argv[i], "-ff") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -ff\n") ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option -ff too small (should be a number/bitmask between 0-255)\n") ;
				usage();
				exit(1);
			}
			OUTPUT_FORMAT_FLAGS = tmp;
			
		}

		//leftover file
		if (strcmp(argv[i], "-u") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -u\n") ;
				usage();
				exit(1);
			}
			i++;
			LEFTOVER_FILE_NAME.assign(argv[i]);
		}

		//verbose
		if (strcmp(argv[i], "-v") == 0) {
			not_defined = 0;
			VERBOSE = 3;
		}

		//scores out
		if (strcmp(argv[i], "-e") == 0) {
			not_defined = 0;
			SCORES_OUT = 0;
		}

		// do not align using reverse index
		if (strcmp(argv[i], "-r") == 0) {
			not_defined = 0;
			MAP_REVERSE = 0;
		}


		if (_personality == Palmapper && strcmp(argv[i], "-S") == 0) {
			not_defined = 0;
			SPLICED_HITS = 1;
		}

		//non-consensus-search
		if (strcmp(argv[i], "-non-consensus-search") == 0) {
			not_defined = 0;
			non_consensus_search = 1;
		}

		// extend the seed-length if too many seed-matches were found
		char const *indexExtentOpt = _personality == Palmapper ? "-index-extend" : "-x";
		if (strcmp(argv[i], indexExtentOpt) == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option %s\n", indexExtentOpt) ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option %s too small\n", indexExtentOpt) ;
				usage();
				exit(1);
			}
			INDEX_DEPTH_EXTRA = tmp ;
		}

		// extend the seed-length if too many seed-matches were found
		char const *indexExtendThresholdOpt = _personality == Palmapper ? "-index-extend-threshold" : "-y";
		if (strcmp(argv[i], indexExtendThresholdOpt) == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option %s\n", indexExtendThresholdOpt) ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option %s too small\n", indexExtendThresholdOpt) ;
				usage();
				exit(1);
			}
			INDEX_DEPTH_EXTRA_THRESHOLD = tmp ;
		}

		// cancel seed processing if more seed matches exist than a threshold
		char const *seedHitCancelThresholdOpt = _personality == Palmapper ? "-seed-hit-cancel-threshold" : "-s";
		if (strcmp(argv[i], seedHitCancelThresholdOpt) == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option %s\n", seedHitCancelThresholdOpt) ;
				usage();
				exit(1);
			}
			i++;
			int tmp = atoi(argv[i]);
			if (tmp < 0) {
				fprintf(stderr, "ERROR: Argument for option %s too small\n", seedHitCancelThresholdOpt) ;
				usage();
				exit(1);
			}
			SEED_HIT_CANCEL_THRESHOLD = tmp ;
		}

		//report all hits-strategy
		if (strcmp(argv[i], "-a") == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			OUTPUT_FILTER=OUTPUT_FILTER_ALL ;
			ALL_HIT_STRATEGY = 1 ;
		}

		char const *outputFilteToprOpt = _personality == Palmapper ? "-z" : "-n";
		if (strcmp(argv[i], outputFilteToprOpt) == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option %s\n", outputFilteToprOpt) ;
				usage();
				exit(1);
			}
			i++;
			
			if ((OUTPUT_FILTER_NUM_TOP = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0' || OUTPUT_FILTER_NUM_TOP<0) 
				{
					fprintf(stderr,
							"ERROR: Number of alignments must be a positive integer value!\n");
					exit(1);
				}
			}
			OUTPUT_FILTER=OUTPUT_FILTER_TOP ;
			SUMMARY_HIT_STRATEGY=1 ;
		}

		//max nr of alignments per read
		/*if (strcmp(argv[i], "-ar") == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -ar\n") ;
				usage();
				exit(1);
			}
			i++;
			if (OUTPUT_FILTER_NUM_LIMIT != 0) {
				fprintf(stderr,
						"ERROR: options -a and -ar exclude themselves!\n");
				exit(1);
			}
			if ((OUTPUT_FILTER_NUM_LIMIT = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of alignments must be an integer value!\n");
					exit(1);
				}
			}
			OUTPUT_FILTER = OUTPUT_FILTER_LIMIT ;
		}*/

		//max nr of alignments per read, randomly chosen!
		char const *outputFilterNumLimitOpt = _personality == Palmapper ? "-n" : "-b";
		if (strcmp(argv[i], outputFilterNumLimitOpt) == 0) {
			not_defined = 0;
			if (OUTPUT_FILTER!=OUTPUT_FILTER_DEFAULT)
			{
				fprintf(stderr, "ERROR: output filter already defined\n") ;
				exit(1) ;
			}
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option %s\n", outputFilterNumLimitOpt) ;
				usage();
				exit(1);
			}
			i++;
			if ((OUTPUT_FILTER_NUM_LIMIT = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of alignments must be an integer value!\n");
					exit(1);
				}
			}
			//OUTPUT_FILTER_NUM_LIMIT = -OUTPUT_FILTER_NUM_LIMIT ;
			OUTPUT_FILTER = OUTPUT_FILTER_LIMIT ;
			BEST_HIT_STRATEGY = 1 ;
		}

		//max number of allowed edit operations
		if (strcmp(argv[i], "-E") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -E\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((NUM_EDIT_OPS = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of edit operations must be an integer value!\n");
					exit(1);
				}
			}
			if (NUM_EDIT_OPS > MAX_EDIT_OPS) {
				fprintf(
						stderr,
						"ERROR: Number of allowed mismatches exceeds maximal number of edit operations (=%d)! Please restart with a lower value!\n",
						MAX_EDIT_OPS);
				exit(1);
			}
		}

		//max number of allowed mismatches
		if (strcmp(argv[i], "-M") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -M\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((NUM_MISMATCHES = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of mismatches must be an integer value!\n");
					exit(1);
				}
			}
			if (NUM_MISMATCHES > MAX_EDIT_OPS) {
				fprintf(
						stderr,
						"ERROR: Number of allowed mismatches exceeds maximal number of edit operations (=%d)! Please restart with a lower value!\n",
						MAX_EDIT_OPS);
				exit(1);
			}
		}

		//max number of allowed gaps
		if (strcmp(argv[i], "-G") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -G\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((NUM_GAPS = atoi(argv[i])) == 0) {
				if (argv[i][0] != '0') {
					fprintf(stderr,
							"ERROR: Number of gaps must be an integer value!\n");
					exit(1);
				}
			}
			if (NUM_GAPS > MAX_EDIT_OPS) {
				fprintf(
						stderr,
						"ERROR: Number of allowed gaps exceeds maximal number of edit operations (=%d)! Please restart with a lower value!\n",
						MAX_EDIT_OPS);
				exit(1);
			}
		}

		//match score
		if (strcmp(argv[i], "-match-score") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -match_score\n") ;
				usage();
				exit(1);
			}
			i++;
			M_SCORE = atof(argv[i]);
		}

		//mismatch score
		if (strcmp(argv[i], "-m") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -m\n") ;
				usage();
				exit(1);
			}
			i++;
			MM_SCORE = atof(argv[i]);
			if (MM_SCORE < 0) {
				fprintf(stderr, "ERROR: Mismatch score must be positive!\n");
				exit(1);
			}
			if (MM_SCORE == 0)
				fprintf(stderr,
						"\n!!! WARNING: Mismatch score is 0! This could lead to bad alignments!\n\n");
		}

		//gap score
		if (strcmp(argv[i], "-g") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -g\n") ;
				usage();
				exit(1);
			}
			i++;
			GAP_SCORE = atof(argv[i]);
			if (GAP_SCORE < 0) {
				fprintf(stderr, "ERROR: Gap score must be positive!\n");
				exit(1);
			}
			if (GAP_SCORE == 0)
				fprintf(stderr,
						"\n!!! WARNING: Gap score is 0! This could lead to bad alignments!\n\n");
		}

		//hitlength limit
		if (strcmp(argv[i], "-l") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -l\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((HITLEN_LIMIT = atoi(argv[i])) == 0 || HITLEN_LIMIT < 0) {
				fprintf(stderr,
						"ERROR: Hitlength limit must be a positive integer value unequal to 0!\n");
				exit(1);
			}
		}

		//min spliced segment length
		if (strcmp(argv[i], "-min-spliced-segment-len") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) 
			  {
			    fprintf(stderr, "ERROR: Argument missing for option -min-spliced-segment-len\n") ;
			    usage();
			    exit(1);
			  }
			i++;
			if ((SPLICED_MIN_SEGMENT_LENGTH = atoi(argv[i])) == 0) {
			  fprintf(stderr,
				  "ERROR: min-spliced-segment-length must be a positive integer value !\n");
			  exit(1);
			}
		}

		//chromosome container size
		if (strcmp(argv[i], "-c") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -c\n") ;
				usage();
				exit(1);
			}
			i++;
			if ((CHROM_CONTAINER_SIZE += atoi(argv[i])) == 0) {
				fprintf(stderr,
						"ERROR: Chromosome Container Size must be an integer value unequal to 0!\n");
				exit(1);
			}
		}

		//last read to be mapped
		if (strcmp(argv[i], "-to") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -to\n") ;
				usage();
				exit(1);
			}
			i++;
			
			if ((LAST_READ_NR = atoi(argv[i])) == 0 || LAST_READ_NR < 0) {
				fprintf(stderr, "ERROR: last read number (-to argument) must be a positive integer value!\n");
				exit(1);
			}
			if (LAST_READ_NR < FIRST_READ_NR) {
				fprintf(stderr, "ERROR: last read number must be greater than first read number (-to and -from options)\n");
				exit(1);
			}
		}

		//first read to be mapped
		if (strcmp(argv[i], "-from") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option -from\n") ;
				usage();
				exit(1);
			}
			i++;
			
			if ((FIRST_READ_NR = atoi(argv[i])) == 0 || FIRST_READ_NR < 0) {
				fprintf(stderr, "ERROR: first read number (-from argument) must be a positive integer value!\n");
				exit(1);
			}
			if (LAST_READ_NR < FIRST_READ_NR) {
				fprintf(stderr, "ERROR: last read number must be greater than first read number (-to and -from options)\n");
				exit(1);
			}
		}

		//nr of unallowed chars in reads
		/*if(strcmp(argv[i],"-n")==0){
		 not_defined = 0;
		 if (i+1 > argc - 1){ usage(); exit(1); }
		 i++;
		 if ((NUM_UNALLOWED_CHARS = atoi(argv[i])) == 0 && argv[i][0] != '0') {
		 fprintf(stderr, "ERROR: Number of non-base symbols in read must be an integer value!\n");
		 exit(1);
		 }
		 }*/

		//gaps most right
		if (strcmp(argv[i], "-d") == 0) {
			not_defined = 0;
			GAPS_MOST_RIGHT = 1;
		}

		//overhang alignment
		if (strcmp(argv[i], "-h") == 0) {
			not_defined = 0;
			OVERHANG_ALIGNMENT = 1;
		}

		//overhang alignment
		if (strcmp(argv[i], "-w") == 0) {
			not_defined = 0;
			STRINGENT_GAPLIMIT = 0;
		}

		//statistics
		if (strcmp(argv[i], "-stat") == 0) {
			not_defined = 0;
			STATISTICS = 1;
		}

		//print out gene, too (for every hit) used for WMD2
		if (strcmp(argv[i], "-target_info") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr,
						"-target_info needs an integer value (1: length only, 2: length+sequence)\n");
				exit(1);
			}
			i++;
			if (((PRINT_SEQ = atoi(argv[i])) == 0) || PRINT_SEQ < 1
					|| PRINT_SEQ > 2) {
				fprintf(stderr, "-target_info value must be either 1 or 2!\n");
				exit(1);
			}
		}

		//flanking region of a hit
		if (strcmp(argv[i], "-flanking") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr,
						"-flanking needs an integer value, the length of one of the symmetric flanking regions of a hit\n");
				exit(1);
			}
			i++;
			if (((FLANKING = atoi(argv[i])) == 0) || FLANKING < 0 || FLANKING
					> 100) {
				fprintf(stderr,
						"-flanking value must be a positive integer and must not exceed 100!\n");
				exit(1);
			}
		}

		//bisulfite mode
		if (strcmp(argv[i], "-B") == 0) {
			not_defined = 0;
			BSSEQ = 1;
		}

		//read2 file
		if (strcmp(argv[i], "-q1") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option %s\n", argv[i]);
				usage();
				exit(1);
			}
			i++;
			Q1_QUERY_FILE_NAMES=strdup(argv[i]);
		}

		if (strcmp(argv[i], "-q2") == 0) {
			not_defined = 0;
			if (i + 1 > argc - 1) {
				fprintf(stderr, "ERROR: Argument missing for option %s\n", argv[i]);
				usage();
				exit(1);
			}
			i++;
			Q2_QUERY_FILE_NAMES=strdup(argv[i]);
		}

		if (_personality == Palmapper) {
			//print out gene, too (for every hit) used for WMD2
			if (strcmp(argv[i], "-qpalma") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr,
							"-qpalma needs an argument\n");
					exit(1);
				}
				i++;
				QPALMA_FILE.assign(argv[i]);
			}

			if (strcmp(argv[i], "-no-qpalma") == 0) {
				not_defined = 0;
				NO_QPALMA = true;
				
			}

			if (strcmp(argv[i], "-acc") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr,
							"-acc needs an argument\n");
					exit(1);
				}
				i++;
				ACC_FILES.assign(argv[i]);
			}

			if (strcmp(argv[i], "-don") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr,
							"-don needs an argument\n");
					exit(1);
				}
				i++;
				DON_FILES.assign(argv[i]);
			}

			if (strcmp(argv[i], "-no-ss-pred") == 0) {
				not_defined = 0;
				NO_SPLICE_PREDICTIONS = 1 ;
			}

			if (strcmp(argv[i], "-qpalma-prb-offset-fix") == 0) {
				not_defined = 0;
				QPALMA_PRB_OFFSET_FIX = 1 ;
			}

			if (strcmp(argv[i], "-acc-consensus") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "-acc-consensus needs an argument\n");
					exit(1);
				}
				i++;
				std::string consensus_list = argv[i] ;
				ACC_CONSENSUS.clear() ;
				split_string(consensus_list, ACC_CONSENSUS, ',') ;
				postprocess_consensus_list(ACC_CONSENSUS) ;
			}

			if (strcmp(argv[i], "-don-consensus") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "-don-consensus needs an argument\n");
					exit(1);
				}
				i++;
				std::string consensus_list = argv[i] ;
				DON_CONSENSUS.clear() ;
				split_string(consensus_list, DON_CONSENSUS, ',') ;
				postprocess_consensus_list(DON_CONSENSUS) ;
			}

			// qpalma-triggered reads file                     //  #A#
			if (strcmp(argv[i], "-log-triggered-reads") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1) {
					fprintf(stderr, "ERROR: Argument missing for option -log-triggered-reads\n") ;
					usage();
					exit(1);
				}
				i++;
				LOG_TRIGGERED = 1;
				TRIGGERED_LOG_FILE.assign(argv[i]);
			}                                                   // #A#

			//report output file
			if (strcmp(argv[i], "-stranded") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1 || (strcmp(argv[i + 1], "left") != 0 && strcmp(argv[i + 1], "right") != 0 && strcmp(argv[i + 1], "plus") != 0 && strcmp(argv[i + 1], "minus") != 0)) {
					fprintf(stderr, "ERROR: Argument missing for option -stranded\nMust be [ left | right | plus | minus ]") ;
					usage();
					exit(1);
				}
				i++;
				if ((strcmp(argv[i], "left") == 0 || strcmp(argv[i], "plus") == 0)) {
					STRAND = 1 ; // plus
				}
				else {
					STRAND = 0 ; // minus
				}
			}

			//Used protocol
			if (strcmp(argv[i], "-protocol") == 0) {
				not_defined = 0;
				if (i + 1 > argc - 1 || (strcmp(argv[i + 1], "first") != 0 && strcmp(argv[i + 1], "second") != 0)) {
					fprintf(stderr, "ERROR: Argument missing for option -protocol\nMust be [ first | second ]") ;
					usage();
					exit(1);
				}
				i++;
				if (strcmp(argv[i], "first") == 0) {
					PROTOCOL = PROTOCOL_FIRST ; // first
				}
				else {
					PROTOCOL = PROTOCOL_SECOND ; // second
				}
			}


			if (not_defined == 1) {
				fprintf(stderr, "ERROR: unknown option %s\n", argv[i]) ;
				usage();
				exit(1);
			}
		}
	}

	//Initialize query file vectors (name and strand information)
	postprocess_query_filenames(Q_QUERY_FILE_NAMES, STRAND) ;
	postprocess_query_filenames(Q1_QUERY_FILE_NAMES, 1) ;
	postprocess_query_filenames(Q2_QUERY_FILE_NAMES, 0) ;

	if (has_index == 0 || QUERY_FILE_NAMES.size() == 0 || has_genome == 0) {
		usage();
		exit(1);
	}

	if (BWA_INDEX && REPORT_REPETITIVE_SEEDS)
	{
		fprintf(stderr, "The combination of the two options -bwa and -report-rep-seed is not implemented yet\n") ;
		exit(-1) ;
	}
	
	NOT_MAXIMAL_HITS = SEED_HIT_CANCEL_THRESHOLD || INDEX_DEPTH_EXTRA_THRESHOLD;

	return 0;
}

void Config::VersionHeader() {
	if (getPersonality() == Palmapper) {
		printf("\nPALMapper version %s   (PALMapper is a fusion of GenomeMapper & QPALMA)\n", VERSION);
		printf("written by Korbinian Schneeberger, Joerg Hagmann, Gunnar Raetsch, Geraldine Jean, Fabio De Bona, Stephan Ossowski, and others\n");
		printf("Max Planck Institute for Developmental Biology and Friedrich Miescher Laboratory, Tuebingen, Germany, 2008-2010\n\n");
	} else {
		printf("\nGenomeMapper version %s\n", VERSION);
		printf("written by Korbinian Schneeberger, Joerg Hagmann, Stephan Ossowski, Felix Ott, Gunnar Raetsch and others\n");
		printf("Max Planck Institute for Developmental Biology, Tuebingen, Germany, 2008-2010\n\n");
	}
}

int Config::usage() {
	//VersionHeader() ;

	if (getPersonality() == Palmapper) {
		printf("USAGE: palmapper [options]\n");

		printf("\n");

		//MANDATORY parameters
		printf("mandatory:\n");
		printf(" -i STRING                       reference sequence (fasta file and prefix to index files)\n");
		printf(" -q STRING[,STRING,..,STRING]    query filename (fasta, fastq, or SHORE flat file)\n");
		printf(" -q1 STRING[,STRING,..,STRING]   \"left\" query filename for paired-end reads (fasta, fastq, or SHORE flat file)\n");
		printf(" -q2 STRING[,STRING,..,STRING]   \"right\" query filename for paired-end reads (fasta, fastq, or SHORE flat file)\n");

		//OPTIONAL parameters
		printf("\n\n");
		printf("optional:\n");

		printf(" -stranded STRING        strand specific experiment (left, right, plus, minus)\n");
		printf(" -protocol STRING        protocol used to prepared RNA-seq data (first,second)\n");
		printf("                         examples: RNA ligation is first and dUTP protocol is second strand\n");
		printf(" -f STRING               output format (\"shore\", \"bed\", \"bedx\", \"sam\", \"bam\", \"bamp\" or \"bamn\")[sam]\n");
		printf(" -samtools STRING        explicit samtools path (used for bam output)\n");
		printf(" -ff INT                     bitwise output sam format flag (0x1: read sequence, 0x2: read quality, 0x4: common sam flags, 0x8: extended same flags)[15]\n");
		printf(" -include-unmapped-reads     write directly unmapped reads in sam file\n");
		printf(" -o STRING                   output filename [stdout]\n");
		printf(" -H STRING                   output filename for spliced hits [no output]\n");
		printf(" -u STRING                   output filename for unmapped reads [/dev/null]\n\n");

		printf(" -rlim INT      limit the number of reads for alignment\n");
		printf(" -from INT      skip the first <from> reads from query file\n");
		printf(" -to   INT      map only the first <to> reads from query file\n\n");

		printf(" -a             report all alignments\n");
		printf(" -ar INT        report a limited number of alignments (random subset) [10]\n");
		printf(" -z  INT        report a number of top alignments [5]\n");
		printf(" -n  INT        report a maximal number of best alignments\n\n");

		printf(" -r             disable alignment on reverse strand [enabled]\n");
		printf(" -h             perform alignment of flanking regions of hits first [whole read alignment]\n");
		printf(" -d             align gaps most right (ignored for spliced alignments) [most left]\n");
 		printf(" -w             allow more gaps for best hit (ignored for spliced alignments) [retain gap limit]\n\n");

		printf(" -bwa INT                         use burrows-wheeler index instead of k-mer index (bwa-based) with a given seed length\n") ;
		printf(" -seed-hit-cancel-threshold INT   number of hits of a seed that lead to its ignoration\n");
		printf(" -index-extend-threshold INT      number of hits of a seed that lead to a seed-length extension\n");
		printf(" -index-extend INT                length of seed-length extension\n");
		printf(" -index-precache                  linearly read index file to fill caches\n");
		printf(" -l INT                           minimal considered hit length [seed length]\n");
		printf(" -c INT                           seed container size [15.000.000]\n\n");

		printf(" -threads INT                     maximal number of threads [1] \n");
		printf(" -v                               verbose [silent]\n\n");

		printf(" -rtrim INT              shortens the read until a hit is found or the minimal length is reached\n");
		printf(" -rtrim-step INT         rtrim step size\n");
		printf(" -polytrim INT           trims polyA or polyT ends until a hit is found or the minimal length is reached\n");
		printf(" -fixtrim INT            shortens the read to a fixed length\n");
		printf(" -fixtrimleft INT        Removes the given number of first nucleotides of each read (can be used with -fixtrimright but not -fixtrim)\n");
		printf(" -fixtrimright INT       Removes the given number of last nucleotides of each read  (can be used with -fixtrimleft but not -fixtrim)\n");
		printf(" -adaptertrim INT        trims away known adapter sequences, read is dropped, when shorter than parameter\n\n");
		//printf(" -adaptertrim-log STRING    \n");


		printf(" -M INT                   max number of mismatches [auto]\n");
		printf(" -G INT                   max number of gaps [auto]\n");
		printf(" -E INT                   max edit operations [auto]\n");
		printf(" -m DOUBLE                mismatch penalty [4]\n");
		printf(" -g DOUBLE                gap penalty [5]\n");
		printf(" -match-score DOUBLE      match penalty [0]\n\n");


		printf(" -S                       report spliced alignments (detailed options below)\n");
		printf("spliced hits definitions: (-S required)\n");
		printf(" -qpalma STRING                        file name with qpalma parameters (essential)\n");
		printf(" -qpalma-use-map-max-len INT           limit the map extension up- and downstream to the given length [10.000]\n");
		printf(" -qpalma-prb-offset-fix                automatically fix the quality offset, if necessary \n\n");

		printf(" -acc STRING                           path name to acceptor splice site predictions (essential if -no-ss-pred not provided)\n");
		printf(" -don STRING                           path name to donor splice site predictions (essential if -no-ss-pred not provided)\n");
		printf(" -acc-consensus STRING                 defines consensus sequences for acceptor sites (separated by \",\") [AG]\n");
		printf(" -don-consensus STRING                 defines consensus sequences for donor sites (separated by \",\") [GT,GC]\n");
		printf(" -no-ss-pred                           indicates that no splice site predictions should be used and only scores positions corresponding to consensus sequences for acceptors and donors\n");
		printf(" -non-consensus-search                 switch on spliced alignments with non consensus sequences as plausible splice sites\n");
		printf(" -score-annotated-splice-sites STRING[,STRING,..,STRING]  set score of annotated splice sites from gff3 files to 1\n");
		printf(" -junction-remapping STRING[,STRING,..,STRING]  enables remapping of unmapped or unspliced reads against the junction list provided in gff3 files\n");
		printf(" -junction-remapping-coverage INT      minimum alignment support to take into account a junction\n\n");

		printf(" -filter-splice-sites-top-perc FLOAT   trigger spliced alignments, if read covers top percentile splice site (between 0 and 1) [0.01]\n");
		printf(" -filter-splice-region INT             extension of the read region up- and downstream for triggeringspliced alignments by presence of splice sites [5]\n");
		printf(" -filter-max-edit INT                  trigger spliced alignment, if unspliced alignment has at least this many edit operations [0]\n");
		printf(" -filter-max-mismatches INT            trigger spliced alignment, if unspliced alignment has at least this many mismatches [0]\n");
		printf(" -filter-max-gaps INT                  trigger spliced alignment, if unspliced alignment has at least this many gaps [0]\n");
		printf(" -log-triggered-reads STRING           log file containing the triggered reads\n\n");

		printf(" -C INT                                min combined length [auto]\n");
		printf(" -L INT                                min length of long hit [auto]\n");
		printf(" -K INT                                min length of short hit [auto]\n");
		printf(" -SA INT                               maximum number of spliced alignments per read [10]\n");
		printf(" -NI INT                               maximum number of introns in spliced alignments [auto]\n");
		printf(" -CT INT                               distance to tolerate between hit and existing hit cluster [10]\n");
		printf(" -QMM INT                              number of matches required for identifying a splice site [auto]\n");
		printf(" -I INT                                longest intron length  [auto]\n");
		printf(" -MI INT                               shortest intron length  [30]\n");
		printf(" -min-spliced-segment-len INT          minimal exon length [auto]\n");
		//printf(" -EL INT                               minimal number of nucleotides in a spliced segment [auto]\n\n") ;

		printf(" -report STRING                        file for map reporting\n");
		printf(" -report-ro STRING                     file for map reporting (read only)\n");
		//printf(" -report-reset                         does not load report even if it is available\n");
		printf(" -report-rep-seed                      switch on reporting of repetitive seeds\n");
		printf(" -report-map-region                    switch on reporting of mapped regions\n");
		//printf(" -no-report-map-region                 switch off reporting of mapped regions\n");
		printf(" -report-map-read                      switch on reporting of mapped reads\n");
		//printf(" -no-report-map-read                   switch off reporting of mapped reads\n");
		printf(" -report-spliced-read                  switch on reporting of spliced reads\n");
		//printf(" -no-report-spliced-read               switch off reporting of spliced reads\n");
		printf(" -report-splice-sites FLOAT            report splice sites with confidence not less that threshold\n");
		printf(" -report-splice-sites-top-perc FLOAT   report splice sites with confidence in top percentile (between 0 and 1)\n");
		printf(" -report-gff-init STRING               initialize map with exons from GFF file\n");
		printf(" -report-coverage-map STRING           report genome coverage in map format\n");
		printf(" -report-coverage-wig STRING           report genome coverage in wiggle format\n\n");
		printf(" -report-junctions STRING              report splice site junctions in gff3 format\n\n");
		//printf(" -qpalma-use-map                       use map for qpalma alignments\n");
		//printf(" -e         report edit operations (alignment scores)\n");

	} else {

		printf("USAGE: genomemapper [options] -i <reference> -q <reads>\n");

		printf("\n");
		printf("*Mandatory arguments*\n\n");
		printf("    -i STRING                     reference sequence (fasta file and prefix to index files)\n");
		printf("    -q STRING[,STRING,..,STRING]  query filename (fasta, fastq, or SHORE flat file)\n");
		printf("\n");

		printf("*Optional arguments*\n\n");
		printf("  Output:\n");
		printf("    -f STRING        output format (\"shore\" or \"sam\") [shore]\n");
		printf("    -o STRING        output filename [stdout]\n");
		printf("    -u STRING        output filename for unmapped reads [no output]\n");

		printf("\n  Mapping strategies (max. one selectable):\n");
		printf("    -a               all-hit-strategy: report all valid alignments\n");
		//printf("    -ar INT          report a limited number of alignments [best alignments only]\n");
		printf("    -n INT           n-best-hit-strategy: report the best n alignments\n");
		printf("    -b INT           best-hit-strategy with limit: report b number of alignments\n                     with top score\n");
		printf("    [default:        best-hit-strategy: report all alignments with top score]\n");

		printf("\n  Alignment settings:\n");
		printf("    -M INT           max number of mismatches [3]\n");
		printf("    -G INT           max number of gaps [1]\n");
		printf("    -E INT           max edit operations [3]\n\n");
		printf("    -m DOUBLE        mismatch penalty [4.0]\n");
		printf("    -g DOUBLE        gap penalty [5.0]\n\n");
		printf("    -d               align gaps most right [most left]\n");
		printf("    -e               report edit operations instead of alignment scores\n");

		printf("\n  Heuristic speedup:\n");
		printf("    -s INT           seed hit cancel threshold: number of hits of a seed that\n                     lead to its ignoration\n");
		printf("    -y INT           seed extend threshold: number of hits of a seed that lead\n                     to a seed-length extension\n");
		printf("    -x INT           seed extend: length of seed-length extension\n");
		printf("    -h               perform alignment of flanking regions of hits first\n                     [whole read alignment]\n");
		printf("    -l INT           seed length: do not align hits shorter than <l> [index depth]\n");

		printf("\n  Read trimming:\n");
		printf("    -fixtrim INT     shortens the read to fixed length\n");
		printf("    -rtrim INT       shortens the read until a hit is found or the minimal length\n                     is reached\n");
		printf("    -rtrim-step INT  trimming step size in bp for -rtrim\n");
		printf("    -polytrim INT    trims polyA or polyT ends until a hit is found or the minimal\n                     length is reached\n");
		printf("    -adaptertrim INT trims away known adapter sequences; read is dropped when\n                     shorter than parameter\n");

		printf("\n  Read subsets:\n");
		printf("    -to INT          map only the first <to> reads from query file\n");
		printf("    -from INT        skip the first <from> reads from query file\n");

		printf("\n  Mapping of bisulfite-treated reads:\n");
		printf("    -B               activates bisulfite mapping mode. There are two modi:\n                     - 2conversion-mode for paired end reads (either -q1 and -q2\n                       or -q and SHORE format)\n                     - 4conversion-mode for single reads (-q)\n");
		printf("    -q1              filename of paired end reads 1 for 2conversion mode\n");
		printf("    -q2              filename of paired end reads 2 for 2conversion mode\n");

		printf("\n  Other options:\n");
		printf("    -r               align only on forward strand of reference [both]\n");
		printf("    -t INT           maximal number of threads [4] \n");
		//printf("    -w               allow more gaps for best hit [retain gap limit]\n");
		printf("    -c INT           seed container size (value will be added to 15Mio) [15.000.000]\n\n");
		printf("    -v               verbose [silent]\n\n");
	}

	return 0;
}
