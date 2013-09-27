#pragma once

#include <palmapper/Genome.h>
#include <palmapper/Hits.h>
#include <palmapper/QPalma.h>
#include <palmapper/Read.h>
#include <palmapper/TopAlignments.h>
#include <palmapper/JunctionMap.h>

class Mapper {

	friend class Hits;

public:

	enum ResultState {
		ReadMapped = 0,
		NothingFound,
		TooShortAfterTrimming,
		MappingFailed,
		NonACGTChar,
		ReadShorterThanHitLengthLimit,
		IgnoreResultBound,
	};
	class Result {
	public:
		Result(Mapper &mapper)
		: 	_orig(mapper._queryFile),
		  	_work(mapper._queryFile),
		  	_readMappings(mapper._genome, mapper._genomeMaps, mapper, _work),
		  	_qpalma(_work, mapper._qpalma)
		{
			_state = NothingFound;
			_rtrim_cut=0;
			_polytrim_cut_start=0;
			_polytrim_cut_end=0;
		}

		ResultState _state;
		Read _orig;
		Read _work;
		Hits _readMappings;
		QPalma::Result _qpalma;

		int _rtrim_cut;
		int _polytrim_cut_start;
		int _polytrim_cut_end;
	};

	class Reporter {
	public:
		virtual void report(Result &result, JunctionMap &junctionmap) = 0;
	};

	Mapper(Genome const &genome, GenomeMaps &genomemaps, QueryFile &queryFile, QPalma &qpalma, Reporter &reporter, JunctionMap &junctionmap, JunctionMap &annotatedjunctions);
	~Mapper();
	void setProgressChar(char c) {
		_progressChar = c;
	}

	int map_reads() ;
	int REDUNDANT;

protected:
	void map_reads_timing(int count_reads, float this_read=-1);
	int init_constants()  ;
	int init_statistic_vars() ;
	int init_operators() ;
	int init_alignment_structures(Config * config);
	void map_read(Result &result, clock_t start_time);
	GenomeArr GENOME;

	Genome const &_genome;
	GenomeMaps &_genomeMaps ;
	JunctionMap & _junctionmap;
	JunctionMap & _annotatedjunctions;

private:
	QueryFile &_queryFile;
	QPalma const &_qpalma;
	Reporter &_reporter;

	int num_spliced_alignments_triggered;
	char _progressChar;

	FILE *_ADAPTERTRIM_LOG_FP;
	FILE *_TRIGGERED_LOG_FP;

	clock_t time1, time2a, time2b, time2c, time3;
	clock_t last_timing_report;

	unsigned int MAX_USED_SLOTS;
	unsigned int MAXHITS;
	int c_map_fast;
	int c_map_short_read;
	std::vector<bool> seed_covered ;
	CHROMOSOME_ENTRY_CONTAINER CHROMOSOME_ENTRY_OPERATOR;
};
