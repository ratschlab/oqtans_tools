#pragma once

#include <palmapper/Chromosome.h>
#include <palmapper/QueryFile.h>

extern int compare_editops(const void *a, const void *b) ;
extern char get_compl_base(char c);
extern void print_stats(QueryFile &queryFile);
extern int print_perfect_hits(unsigned int num);
extern int print_largest_hit();
extern void print_leftovers(Read const &read, const char *tag, std::ostream *LEFTOVER_FP);
extern void print_alignment_matrix(Read const &read, int chrstart, int readstart, int length, int offset_front, int offset_end, Chromosome const &chr, char ori, int K);
extern int compare_int(const void *a, const void *b) ;
