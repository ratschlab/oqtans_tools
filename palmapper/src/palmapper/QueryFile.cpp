#include <iostream>

#include <palmapper/QueryFile.h>
#include <palmapper/Read.h>

using namespace std;

QueryFile::QueryFile(std::vector<std::string> const &filenames,std::vector<int> const &strands) : _filenames(filenames), _strands(strands) {
	_currentFile = -1;
	_current_strand=-1;
	_readCount = 0;
	_file = NULL;
	open_next_file();
}


QueryFile::~QueryFile() {
	::fclose(_file);
}

bool QueryFile::next_line(char *buf, int maxLen) {
	if (::fgets(buf, maxLen, _file) == NULL)
		return false;
	++_lineNr;
	return true;
}

Read *QueryFile::next_read() {
	Read *read = new Read(*this);
	if (!next_read(*read)) {
		delete read;
		return NULL;
	}
	return read;
}

bool QueryFile::next_read(Read &read, int &strand) {
	Mutex::Locker locker(_mutex);
	while (read.read_short_read() > 0) {
		if (!open_next_file()) {
			if (_readCount == 0)
				cerr << "\n!!! WARNING: None of the given file(s) contain any usable read!\n\n";
			return false;
		}
	}
	strand=_current_strand;
	read._nr = _readCount++;
	if (read.length() > _maxReadLen)
		_maxReadLen = read.length();
	return true;
}

bool QueryFile::next_read(Read &read) {
	Mutex::Locker locker(_mutex);
	while (read.read_short_read() > 0) {
		if (!open_next_file()) {
			if (_readCount == 0)
				cerr << "\n!!! WARNING: None of the given file(s) contain any usable read!\n\n";
			return false;
		}
	}
	read._nr = _readCount++;
	if (read.length() > _maxReadLen)
		_maxReadLen = read.length();
	return true;
}


bool QueryFile::open_next_file() {
	if (++_currentFile >= _filenames.size())
		return false;
	if (_file != NULL)
		::fclose(_file);
	_file = Util::openFile(_filenames[_currentFile], "r");
	_lineNr = 0;
	_maxReadLen = 0;
	_current_strand=_strands[_currentFile];
	return true;
}

int QueryFile::determine_read_length(std::vector<std::string> const &filenames,std::vector<int> const &strands) {
	QueryFile file(filenames,strands);
	int const sample_size = 10000;
	int sum_read_length = 0;
	int nr_read = 0;

	Read r(file);
	for (;  nr_read < sample_size; ++nr_read) {
		if (!file.next_read(r))
			break;
		sum_read_length += r.length();
	}
	return sum_read_length / nr_read;
}

void QueryFile::Location::printOn(std::ostream &out) const {
	out << *_filename;
	out << '(';
	out << _lineNr;
	out << ')';
}
