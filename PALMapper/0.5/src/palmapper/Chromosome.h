#pragma once

#include <ctype.h>
#include <stdint.h>
#include <string.h>

#include <stdio.h>
#include <string>

template <int bits> class ChromosomeBase {
	friend class Genome;
public:
	ChromosomeBase() : _desc("<undefined>"){
		_length = 0;
		_data = NULL;
		_nr = -1;
	}
	ChromosomeBase(char *ascii, uint32_t len) {
		data(ascii, len);
	}
	~ChromosomeBase() {
		delete[] _data;
	}

	uint32_t length() const {
		return _length;
	}

	uint32_t nr() const {
		return _nr;
	}

	char operator[](uint32_t index) const {
		uint8_t compound = _data[index / _valuesPerByte];
		return num2dna((compound >> ((index % _valuesPerByte) * bits)) & _mask);
	}

	char const *desc() const {
		return _desc.c_str();
	}

	void desc(char const *d) {
		_desc.assign(d);
	}

	void data(char *ascii, uint32_t len);

private:
	static char num2dna(uint8_t num) {
		return _num2dna[num];
	}

	static bool _classInitialized;
	static bool initClass();

	static uint8_t const _mask = (uint8_t) '\xFF' >> (8 - bits);
	static int const _valuesPerByte = 8 / bits;
	static char const *_num2dna;
	static uint8_t _dna2num[256];

	uint32_t _length;
	uint32_t _nr;
	std::string _desc;
	uint8_t *_data;
};

template <int bits> uint8_t ChromosomeBase<bits>::_dna2num[256];

template <> inline char ChromosomeBase<8>::num2dna(uint8_t num) {
	return (char) num;
}

template <> inline void ChromosomeBase<8>::data(char *ascii, uint32_t len) {
	_length = len;
	_data = (uint8_t*) ascii;
}

template <int bits> bool ChromosomeBase<bits>::_classInitialized = ChromosomeBase<bits>::initClass();
template <int bits> bool ChromosomeBase<bits>::initClass() {
	int nrBases = ::strlen(_num2dna);
	for (int i = 0; i < nrBases; ++i) {
		_dna2num[(int)_num2dna[i]] = i;
		_dna2num[::tolower(_num2dna[i])] = i;
	}
	return true;
}

template <int bits> void ChromosomeBase<bits>::data(char *ascii, uint32_t len) {
	if (_classInitialized)
		initClass();
	_length = len;
	_data = new uint8_t[(len + _valuesPerByte - 1) / _valuesPerByte];
	if (len == 0) {
		delete[] ascii;
		return;
	}

	char *src = ascii;
	char *end = ascii + len;
	uint8_t *pos = _data;

	for (;;) {
		uint8_t b = 0;
		for (int i = 0; i < _valuesPerByte; ++i) {
			b |= _dna2num[(unsigned)*src]  << (bits * i);
			if (++src >= end) {
				*pos = b;
				delete[] ascii;
				return;
			}
		}
		*pos++ = b;
	}
}

#ifndef GM
typedef ChromosomeBase<2> Chromosome;
#else
typedef ChromosomeBase<4> Chromosome;
#endif
