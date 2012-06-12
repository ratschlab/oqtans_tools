// authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

// Authors: Gunnar R\"atsch and Lisa Thalheim
// Copyright (C) 2008 by Friedrich Miescher Laboratory of the Max Planck Society

#include "palmapper.h"
#include "print.h"

#include <lang/Thread.h>
#include <palmapper/FileReporter.h>
#include <palmapper/Mapper.h>
#include <palmapper/JunctionMap.h>
#include <wait.h>

Config _config;
Statistics _stats;

using namespace lang;

 
class MapperThread : public Thread, public Mapper {
public:
MapperThread(Genome &genome,	GenomeMaps &genomemaps, QueryFile &queryFile, QPalma &qpalma, Reporter &reporter, JunctionMap &junctionmap, JunctionMap &annotated_junctions)
	: 	Mapper(genome, genomemaps, queryFile, qpalma, reporter,junctionmap,annotated_junctions) {}
	void run() {
		map_reads();
	}
};


int main(int argc, char *argv[]) 
{
	_config.VersionHeader() ;

	init_shogun() ;

	// timing //////////////
	time_t timer, timer_mid, timer2;
  	timer=time(NULL);
  	//printf("The current time is %s",asctime(localtime(&timer)));
	////////////////////////

	// initialize variables
	_config.parseCommandLine(argc, argv);
	

	Genome genome;

	_config.applyDefaults(&genome) ;
	_config.checkConfig() ;

	FILE *OUT_FP = NULL ;
	if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
		OUT_FP = Util::openFile(_config.OUT_FILE_NAME, "w");
	else
	{
		OUT_FP = TopAlignments::open_bam_pipe(_config.OUT_FILE_NAME) ;
		
		if (OUT_FP==NULL)
		{
			fprintf(stderr, "popen failed: samtools subprocess could not be started (mapped reads)\n") ;
			exit(-1) ;
		}
	}
	
	FILE *SP_OUT_FP = OUT_FP ;
	if (_config.SPLICED_HITS && _config.SPLICED_OUT_FILE_NAME.length() > 0 && _config.SPLICED_OUT_FILE_NAME!=_config.OUT_FILE_NAME )
	{
		if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
			SP_OUT_FP = Util::openFile(_config.SPLICED_OUT_FILE_NAME, "w") ;
		else
		{
			SP_OUT_FP = TopAlignments::open_bam_pipe(_config.SPLICED_OUT_FILE_NAME) ;
			
			if (SP_OUT_FP==NULL)
			{
				fprintf(stderr, "popen failed: samtools subprocess could not be started (spliced reads)\n") ;
				exit(-1) ;
			}
		}
	}
    FILE *LEFTOVER_FP = _config.LEFTOVER_FILE_NAME.length() > 0 ? Util::openFile(_config.LEFTOVER_FILE_NAME, "w+") : NULL;


	QueryFile queryFile(_config.QUERY_FILE_NAMES,_config.QUERY_FILE_STRANDS);
	FileReporter reporter(OUT_FP, SP_OUT_FP, LEFTOVER_FP);
	JunctionMap junctionmap(genome,_config.MAP_JUNCTIONS_COVERAGE,_config.ACC_CONSENSUS,_config.DON_CONSENSUS,_config.ACC_CONSENSUS_REV,_config.DON_CONSENSUS_REV);
	JunctionMap annotated_junctions(genome,_config.MAP_JUNCTIONS_COVERAGE,_config.ACC_CONSENSUS,_config.DON_CONSENSUS,_config.ACC_CONSENSUS_REV,_config.DON_CONSENSUS_REV);	
	
	if ( _config.OUTPUT_FORMAT == OUTPUT_FORMAT_BAM)
	{
		TopAlignments::print_bam_header(genome, OUT_FP) ;
		if (SP_OUT_FP!=OUT_FP)
			TopAlignments::print_bam_header(genome, SP_OUT_FP) ;
	}

	if (_config.MAP_JUNCTIONS){
		int ret=junctionmap.init_from_gffs(_config.MAP_JUNCTIONS_FILE);
		if (ret!=0)
			return -1;
	}
	
	if (_config.SCORE_ANNOTATED_SPLICE_SITES){
		int ret=annotated_junctions.init_from_gffs(_config.ANNOTATED_SPLICE_SITES_FILE);
		if (ret!=0)
			return -1;
	}

	// initialize GenomeMaps
//	if (_config.REPORT_REPETITIVE_SEEDS || _config.REPORT_MAPPED_REGIONS || _config.REPORT_MAPPED_READS || _config.REPORT_FILE!=NULL || _config.FILTER_BY_SPLICE_SITES || _config.QPALMA_USE_SPLICE_SITES)
//	{
//		genomemaps.init_reporting() ;
//
//		if (!_config.REPORT_RESET)
//		{
//			genomemaps.read_reporting() ;
//			genomemaps.do_reporting(1) ;
//		}
//	}

	GenomeMaps *genomemaps = NULL;
	QPalma *qpalma = NULL;
	genomemaps = new GenomeMaps(genome);
	qpalma = new QPalma(&genome, genomemaps, 0);

	
	if (_config.SPLICED_HITS && _config.FILTER_BY_SPLICE_SITES && !_config.NO_SPLICE_PREDICTIONS)
	{
		//genomemaps = new GenomeMaps(genome);
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "Using splice sites with confidence in top %1.2f%% percentile for read filtering\n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC) ;
		else
			fprintf(stdout, "Using splice sites with confidence => %1.2f%%  for read filtering\n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_ACC) ;

		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.FILTER_BY_SPLICE_SITES_THRESH_ACC = _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC ;
		if (_config.ACC_FILES.length()>0)
			qpalma->map_splice_sites(_config.ACC_FILES, 'a', _config.FILTER_BY_SPLICE_SITES_THRESH_ACC, _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0, true) ;
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, " -> acceptor  splice sites with confidence >= %1.2f%% \n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_ACC) ;

		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.FILTER_BY_SPLICE_SITES_THRESH_DON = _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC ;
		if (_config.DON_FILES.length()>0)
			qpalma->map_splice_sites(_config.DON_FILES, 'd', _config.FILTER_BY_SPLICE_SITES_THRESH_DON, _config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0, true) ;
		if (_config.FILTER_BY_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, " -> donor  splice sites with confidence >= %1.2f%% \n", 100*_config.FILTER_BY_SPLICE_SITES_THRESH_DON) ;
	}

	if (_config.SPLICED_HITS && _config.QPALMA_USE_SPLICE_SITES && !_config.NO_SPLICE_PREDICTIONS)
	{
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "Using splice sites with confidence in top %1.2f%% percentile to define QPALMA regions\n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC) ;
		else
			fprintf(stdout, "Using splice sites with confidence >= %1.2f to define QPALMA regions\n", _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC) ;
		
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.QPALMA_USE_SPLICE_SITES_THRESH_ACC = _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC ;
		if (_config.ACC_FILES.length()>0)
			qpalma->map_splice_sites(_config.ACC_FILES, 'a', _config.QPALMA_USE_SPLICE_SITES_THRESH_ACC, _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0, false) ;
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "-> acceptor splice sites build_index with confidence >= %1.2f%% \n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_ACC) ;
			
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			_config.QPALMA_USE_SPLICE_SITES_THRESH_DON = _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC ;
		if (_config.DON_FILES.length()>0)
			qpalma->map_splice_sites(_config.DON_FILES, 'd', _config.QPALMA_USE_SPLICE_SITES_THRESH_DON, _config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0, false) ;
		if (_config.QPALMA_USE_SPLICE_SITES_THRESH_TOP_PERC!=0.0)
			fprintf(stdout, "-> donor splice sites with confidence >= %1.2f%% \n", 100*_config.QPALMA_USE_SPLICE_SITES_THRESH_DON) ;
	}

	// slightly hacky ... 
	if (qpalma)
		Read::set_quality_offset(qpalma->get_qpalma_quality_offset()) ;
	
//	if (_config.REPORT_GFF_FILE_NAME.size()>0)
//		genomemaps.init_with_gff(_config.REPORT_GFF_FILE_NAME) ;
	
	// timing //////////////
	timer_mid=time(NULL);
	//printf("The current time is %s",asctime(localtime(&timer_mid)));
  	////////////////////////

 	if (_config.VERBOSE) { printf("Mapping reads\n"); }

 	unsigned int numThreads = _config.NUM_THREADS;
	MapperThread *threads[numThreads];
	std::string threadIds(".+-:=!$'");
	for (unsigned int i = 0; i < numThreads; ++i) {
		threads[i] = new MapperThread(genome, *genomemaps, queryFile, *qpalma, reporter, junctionmap, annotated_junctions);
		threads[i]->setProgressChar(threadIds[i % threadIds.length()]);
		printf("Starting thread %d\n", i);
		threads[i]->launch();
		//threads[i]->run();
	}
	for (unsigned int i = 0; i < numThreads; ++i) {
	        threads[i]->join();
		delete threads[i];
	}
	reporter.done();
	
	if (_config.OUT_FILE_NAME.length() > 0)
	{
		if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
			fclose(OUT_FP);
		else
		{
			int ret=TopAlignments::close_bam_pipe(OUT_FP) ;
			if (ret!=0)
			{
				fprintf(stderr, "samtools pipe failed (ret=%i)\n", ret) ;
				exit(-1) ;
			}
		}
	}
	if (_config.SPLICED_OUT_FILE_NAME.length() > 0 && SP_OUT_FP!=OUT_FP)
	{
		if ( _config.OUTPUT_FORMAT != OUTPUT_FORMAT_BAM )
			fclose(SP_OUT_FP);
		else
		{
			int ret=TopAlignments::close_bam_pipe(SP_OUT_FP) ;
			if (ret!=0)
			{
				fprintf(stderr, "samtools pipe failed (ret=%i)\n", ret) ;
				exit(-1) ;
			}
		}
	}
	if (_config.LEFTOVER_FILE_NAME.length() > 0) 
		fclose(LEFTOVER_FP);

	if (_config.STATISTICS)	{
		print_stats(queryFile);
//TODO:		printf("R E D U N D A N T : %d\n", mapper.REDUNDANT);
//		printf("Max used slots: %d\n", MAX_USED_SLOTS);
		printf("\nList iterations: %d\n", _stats.listcount);
		printf("List iterations occurrences: %d\n",_stats.listocc);
	}
	else if (_config.VERBOSE) printf("Mapped Reads: %i of all %d reads\n", _stats.READS_MAPPED, _stats.NUM_READS);

	// timing //////////////
  	timer2=time(NULL);
  	//printf("The current time is %s",asctime(localtime(&timer2)));
  	int seconds = difftime(timer2, timer_mid);
  	int hours = seconds / 3600;
  	seconds %= 3600;
  	int minutes = seconds / 60;
  	seconds %= 60;
  	if (_config.STATISTICS) {
  		printf("Time needed to pre-process: %gs\n",difftime(timer_mid, timer));
  		printf("Time needed to map: %dh %dm %ds\n", hours, minutes, seconds);
  	}
  	seconds = difftime(timer2, timer);
  	hours = seconds / 3600;
  	seconds %= 3600;
  	minutes = seconds / 60;
  	seconds %= 60;
  	if (_config.STATISTICS) printf("Total time needed: %dh %dm %ds\n", hours, minutes, seconds);
  	////////////////////////

	if (genomemaps != NULL && (_config.REPORT_REPETITIVE_SEEDS || _config.REPORT_MAPPED_REGIONS || _config.REPORT_MAPPED_READS || _config.REPORT_FILE!=NULL))
	{
		genomemaps->do_reporting(1) ;
		genomemaps->write_reporting() ;
		//genomemaps.clean_reporting() ;
	}
	if (_config.REPORT_GENOME_COVERAGE)
		genomemaps->write_cov_reporting() ;

	if (_config.REPORT_JUNCTIONS)
		junctionmap.report_to_gff(_config.REPORT_JUNCTIONS_FILE);
	
	if (qpalma != NULL)
		delete qpalma;
	if (genomemaps != NULL)
		delete genomemaps;

	if (_config.VERBOSE) { printf("Mapping finished\n"); }

	return 0;
}
