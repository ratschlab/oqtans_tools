#pragma once

#include <stdio.h>

#include <vector>

#include <lang/Thread.h>


using lang::Mutex;

class Read;

class QueryFile {

public:

	class Location {
	public:
		Location(std::string const &filename, unsigned long lineNr)
		: _filename(&filename), _lineNr(lineNr) {
		}

		void printOn(std::ostream &out) const;
	private:
		std::string const *_filename;
		unsigned long _lineNr;
	};

	QueryFile(std::vector<std::string> const &filename, std::vector<int> const &strands);
	~QueryFile();

	Read *next_read();

	bool next_read(Read &read, int &strand);
	bool next_read(Read &read);
	bool next_line(char *buf, int maxLen);

	unsigned long line_nr() const {
		return _lineNr;
	}

	Location const getLocation() {
		return Location(_filenames[_currentFile], _lineNr);
	}

	unsigned int max_length() const {
		return _maxReadLen;
	}

	int read_count() const {
		return _readCount;
	}

	void reset_read_count() {
		_readCount=0;
	}

	static int determine_read_length(std::vector<std::string> const &filenames,std::vector<int> const &strands);

private:

	bool open_next_file();

	FILE *_file;
	int _lineNr;
	unsigned int _maxReadLen;
	std::vector<std::string> const &_filenames;
	std::vector<int> const &_strands;

	unsigned long _currentFile;
	int _readCount;
	int _current_strand;
	Mutex _mutex;
};

inline std::ostream &operator<<(std::ostream &out, QueryFile::Location const &loc) {
	loc.printOn(out);
	return out;
}
