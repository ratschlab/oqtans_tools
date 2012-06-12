#include <palmapper/Chromosome.h>

template <>  char const *ChromosomeBase<2>::_num2dna = "ACGT";
template <>  char const *ChromosomeBase<4>::_num2dna = "ACGTNRYMKWSBDHVN";

template <> bool ChromosomeBase<8>::initClass() {
	return true;
}

#if 0
#include <assert.h>
#include <iostream>

extern void func(char a);
extern void func(char a, char b, char c);

int main(int argc, char *argv[]) {
	char const *testData = "ACGTTGCACATG";
	for (unsigned int i = 0; i < ::strlen(testData); ++i) {
		ChromosomeBase<2> c2(::strdup(testData), i);
		ChromosomeBase<4> c4(::strdup(testData), i);
		ChromosomeBase<8> c8(::strdup(testData), i);
		for (unsigned int j = 0; j < i; ++j) {
			assert(c2[j] == testData[j]);
			assert(c4[j] == testData[j]);
			assert(c8[j] == testData[j]);
		}
	}
}

#endif
