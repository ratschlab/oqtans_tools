#include <assert.h>
#include <sstream>

#include <palmapper/FileReporter.h>
#include <palmapper/print.h>

FileReporter::FileReporter(FILE *out, FILE *sp_out, FILE *left_overs) {
	_out = out;
	_sp_out = sp_out;
	_left_overs = left_overs;
	_lastResult = -1;
}

void FileReporter::report(Mapper::Result &result, JunctionMap &junctionmap) {
//	fprintf(stdout,"Delivering result %i\n", result._orig.getNr());

	int readNr = result._work.getNr();
	std::stringstream out;
	std::stringstream sp_out;
	std::stringstream leftOvers;

	if (result._state == Mapper::ReadMapped || (result._state < Mapper::IgnoreResultBound && _config.INCLUDE_UNMAPPED_READS_SAM)){
		result._readMappings.topAlignments().end_top_alignment_record(result._work, &out, &sp_out, result._rtrim_cut, result._polytrim_cut_start, result._polytrim_cut_end,junctionmap);
	} else {
		if (result._state < Mapper::IgnoreResultBound && _config.LEFTOVER_FILE_NAME.length() > 0) {
			char const *text = "";
			switch (result._state) {
				case Mapper::MappingFailed:
					text = " (read mapping failed)";
					break;
				case Mapper::TooShortAfterTrimming:
					text = "(too short after trimming)";
					break;
				default:
					;
			}

			if (!_config.INCLUDE_UNMAPPED_READS_SAM)
				print_leftovers(result._work, text, &leftOvers);
		}
	}
	delete &result;

	Mutex::Locker locker(_mutex);
	while (readNr >= _lastResult + _nrResults) {
		printf("Warning: small result buffer may degrade performance\n");
		_roomLeft.wait(_mutex);
	}
	int index = readNr % _nrResults;
	Entry &e = _results[index];
	assert(!e._used);
	_results[index]._used = true;

	e._out = out.str();
	e._sp_out = sp_out.str();
	e._left_overs = leftOvers.str();

	for (int i = _lastResult + 1;  ; ++i) {
		int pos = i % _nrResults;
		if (!_results[pos]._used)
			break;
		//printf("Writing result %i\n", i);
		Entry &en(_results[pos]);
		assert(en._used);
		print(_out, en._out);
		print(_sp_out, en._sp_out);
		print(_left_overs, en._left_overs);
		en._used = false;
		_lastResult = i;
	}
	_roomLeft.notifyAll();
}


void FileReporter::print(FILE *file, std::string &s){
	if (s.length() == 0 )
		return;
	if (::fwrite(s.c_str(), 1, s.length(), file) != s.length()) {
		fprintf(stderr, "Cound not wrote result\n");
		exit(1);
	}
	s.assign("");
}
