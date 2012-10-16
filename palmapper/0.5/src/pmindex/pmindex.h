// Authors: Korbinian Schneeberger, Stephan Ossowski and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#if defined CONF
#else
#define CONF

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define VERSION "0.3.0"

// ##############################################################
// ####### GLOBAL VARIABLES #####################################
// ##############################################################

#define MIN_INDEX_DEPTH 5
#define MAX_INDEX_DEPTH 15

extern char VERBOSE;
extern char HAS_SLOT;
extern unsigned int SLOT;

extern int INDEX_DEPTH;
extern long int POWER[MAX_INDEX_DEPTH];
extern long int BINARY_CODE[5];

extern unsigned int LONGEST_CHROMOSOME;

extern int debug;

// ##############################################################
// ####### FILE HANDLING ########################################
// ##############################################################

extern char GENOME_FILE_NAME[500];
extern char CHR_INDEX_FILE_NAME[500];
extern char MAPFWD_INDEX_FILE_NAME[500];
extern char META_INDEX_FILE_NAME[500];
extern char GENOME_OUT_FILE_NAME[500];
//char OCC_FILE_NAME[500];

extern FILE *CHR_INDEX_FP;
extern FILE *MAPFWD_INDEX_FP;
extern FILE *META_INDEX_FP;
extern FILE *GENOME_FP;
extern FILE *GENOME_OUT_FP;
//FILE *OCC_FP;

// ##############################################################
// ######## SEQUENCE STORAGE ####################################
// ##############################################################

#define CHR_DESC_LENGTH 50

extern unsigned int CHR_LENGTH;
extern char *CHR_SEQ;

extern char CHR_DESC[2000];
extern char CHR_DESC_TMP[2000];

extern unsigned long int POSITION_COUNTER;

// ##############################################################
// ####### INDEX MANAGEMENT #####################################
// ##############################################################

#define BIN_SIZE 3
#define BIN_SIZE_EXT 20
//#define INDEX_SIZE 16777216 //4^12
//#define INDEX_SIZE 67108864 //4^13
#define INDEX_SIZE 268435456 //4^14
//#define INDEX_SIZE (268435456*4) //4^15
// #define INDEX_SIZE 244140625 // 5^12

#define BLOCK_TABLE_SIZE 16777216	// 2^24 (3 Byte)
#define BLOCK_SIZE 256	// 2^8 (1 Byte)

typedef struct id_structure {
	unsigned char id[4];
} ID;

typedef struct position_structure {
	unsigned int pos;
	unsigned int chr;
} POS;

extern POS *BLOCK_TABLE;
extern unsigned int BLOCK;
extern unsigned short int POSITION;

typedef struct bin_extension_structure {
	ID ids[BIN_SIZE_EXT];
	struct bin_extension_structure *bin_ext;
} BIN_EXT;

typedef struct bin_structure {
	unsigned int num_pos;
	ID ids[BIN_SIZE];
	BIN_EXT *bin_ext;
	BIN_EXT *last_bin_ext;
} BIN;

//BIN *INDEX;
//BIN *INDEX_REV;

//extern BIN *INDEX[INDEX_SIZE];
extern BIN **INDEX;//[INDEX_SIZE];

extern long int NUM_USED_SLOTS; //different to SLOT_COUNTER! This counts the number of different used slots, used by reverse and forward Index.
extern int *USED_SLOTS;//[INDEX_SIZE];


// ##############################################################
// ####### MEMORY MANAGEMENT ####################################
// ##############################################################

#define BIN_EXT_PER_NUGGET 100000

typedef struct nugget_str {
	BIN_EXT buffer[BIN_EXT_PER_NUGGET];
	struct nugget_str *next;
} NUGGET;

typedef struct stmg_str {
	int curr_num;
	NUGGET *nuggets;
} STORAGE;

extern STORAGE *MEM_MGR;

// ##############################################################
// ####### ROUTINES #############################################
// ##############################################################

//init.c
extern int init(int argc, char *argv[]);

//usage.c
extern int usage();

//load.c
extern int load_chromosomes();
extern int desc_parsing(char *c);

//indec.c
extern int index_chromosome(unsigned int chr);

//alloc.c
extern int alloc_bin(int slot);
extern BIN_EXT *alloc_bin_ext() ;
extern void alloc_blocktable();
extern int dealloc_chr();

//write.c
extern int write_meta_index(unsigned int num_chr);
//int write_index();
extern int write_chr_desc(unsigned int chr);

//printindex.c
extern void printindex();

#endif
