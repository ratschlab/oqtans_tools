#pragma once

#include <palmapper/Chromosome.h>
#include <palmapper/Hits.h>
#include <palmapper/Read.h>

#define DIAGONAL 'D'
#define LEFT 'L'
#define UP 'U'

extern Config _config;

extern double WORST_SCORE;
extern double WORST_MM_SCORE;

int kbound_overhang_alignment(Read const &read, HIT* hit, int offset, int readstart, int start, int end, unsigned short int hitreadpos, Chromosome const &chromosome, char orientation, unsigned char mismatches);
int kbound_global_alignment(Read const &read, HIT* hit, unsigned short int hitreadpos, unsigned int start, unsigned int end, Chromosome const &chromosome, char orientation, int Num_edit_ops);


//align.c
int check_mm(Read const &read, Chromosome const &chr, int genome_pos, int readpos, int ori, char conversion);
