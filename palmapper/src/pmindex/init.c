// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "pmindex.h"

int init_defaults();
int init_opts(int argc, char *argv[]);
int init_index();
int init_constants();
int init_genome_file();
int init_chr_index_file();
int init_mapindex_file();
int init_meta_index_file();
int init_genome_out_file();
//int init_stat_file();
int init_mem_master();

int init(int argc, char *argv[]) 
{
	init_defaults();
	init_opts(argc, argv);
	//init_index();
	init_constants();
	init_genome_file();
	init_chr_index_file();
	init_mapindex_file();
	//if (strlen(OCC_FILE_NAME) > 0) init_stat_file();
	if (strlen(GENOME_OUT_FILE_NAME) > 0) init_genome_out_file();
	init_meta_index_file();
	init_mem_master();
	
	alloc_blocktable();

	return 0;
}

int init_defaults() 
{
	int i;

	NUM_USED_SLOTS = 0;
	
	VERBOSE = 0;
	POSITION_COUNTER = 0;
	INDEX_DEPTH = 12;
	LONGEST_CHROMOSOME = 0;

	POWER[0] = 1;
	for(i = 1 ; i < MAX_INDEX_DEPTH; i++) {
		POWER[i] = POWER[i-1] * 4;
	}

	return 0;
}

int init_opts(int argc, char *argv[]) 
{
	int i;
	char has_genome = 0;
	char has_index = 0;

	for (i = 1; i < argc; i++) {

		//genome IN
		if(strcmp(argv[i],"-i")==0){
			if(i+1 > argc - 1){ usage(); exit(1); }
			i++;
			strcpy(GENOME_FILE_NAME, argv[i]);
			has_genome = 1;
			has_index = 1;

			strcpy(CHR_INDEX_FILE_NAME, argv[i]);
                        strcpy(CHR_INDEX_FILE_NAME + strlen(argv[i]), ".cid");

                        strcpy(MAPFWD_INDEX_FILE_NAME, argv[i]);
                        strcpy(MAPFWD_INDEX_FILE_NAME + strlen(argv[i]), ".mfd");

                        strcpy(META_INDEX_FILE_NAME, argv[i]);
                        strcpy(META_INDEX_FILE_NAME + strlen(argv[i]), ".mta");


		}
		
		//genome OUT
		/*if(strcmp(argv[i],"-go")==0){
			if(i+1 > argc - 1){ usage(); exit(1); }
			i++;
			strcpy(GENOME_OUT_FILE_NAME, argv[i]);
			has_genome++;
		}*/
	
/**
		//index
		if(strcmp(argv[i],"-x")==0){
			if(i+1 > argc - 1){ usage(); exit(1); }
			i++;

			strcpy(CHR_INDEX_FILE_NAME, argv[i]);
			strcpy(CHR_INDEX_FILE_NAME + strlen(argv[i]), ".cid");

			strcpy(MAPFWD_INDEX_FILE_NAME, argv[i]);
			strcpy(MAPFWD_INDEX_FILE_NAME + strlen(argv[i]), ".mfd");

			strcpy(MAPREV_INDEX_FILE_NAME, argv[i]);
			strcpy(MAPREV_INDEX_FILE_NAME + strlen(argv[i]), ".mrc");

			strcpy(META_INDEX_FILE_NAME, argv[i]);
			strcpy(META_INDEX_FILE_NAME + strlen(argv[i]), ".mta");

			has_index = 1;
		}
**/

		// seed occurrence stat output
/*		if(strcmp(argv[i],"-o")==0){
			if(i+1 > argc - 1){ usage(); exit(1); }
			i++;
			strcpy(OCC_FILE_NAME, argv[i]);
			has_meta_index = 1;
		}*/

		//depth
		if(strcmp(argv[i],"-s")==0){
			if(i+1 > argc - 1){ usage(); exit(1); }
			i++;
			if ((INDEX_DEPTH = atoi(argv[i])) == 0) {
				fprintf(stderr, "ERROR: seedlength must be an integer value between 5 and 12 (both inclusive)!\n");
				exit(0);
			}
			if (INDEX_DEPTH < MIN_INDEX_DEPTH) {
				fprintf(stderr, "ERROR: specified index depth too low. Minimum allowed is 5!\n");
				exit(1);
			}
			else if (INDEX_DEPTH > MAX_INDEX_DEPTH) {
				fprintf(stderr, "ERROR: specified index depth too big. Maximum allowed is 12!\n");
				exit(1);
			} 
		}
	
		//verbose
		if(strcmp(argv[i],"-v")==0){
			VERBOSE = 1;
		}

  	}

	if (has_index == 0 || has_genome == 0) {
		usage(); exit(1);
	}

	return 1;
}

/*
int init_index()
{
	INDEX_SIZE = 4^INDEX_DEPTH;
	
	INDEX = (BIN *) malloc(INDEX_SIZE * sizeof(BIN));
	INDEX_REV = (BIN *) malloc(INDEX_SIZE * sizeof(BIN));
	
	return 1;
}*/


int init_constants() 
{
	if (INDEX_DEPTH == 5) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 256;		//binary number: 0000 0000 0000 0001 0000 0000
		BINARY_CODE[2] = 512;		//binary number: 0000 0000 0000 0010 0000 0000
		BINARY_CODE[3] = 768;		//binary number: 0000 0000 0000 0011 0000 0000
		BINARY_CODE[4] = 1023;		//binary number: 0000 0000 0000 0011 1111 1111
	}
	if (INDEX_DEPTH == 6) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 1024;		//binary number: 0000 0000 0000 0100 0000 0000
		BINARY_CODE[2] = 2048;		//binary number: 0000 0000 0000 1000 0000 0000
		BINARY_CODE[3] = 3072;		//binary number: 0000 0000 0000 1100 0000 0000
		BINARY_CODE[4] = 4095;		//binary number: 0000 0000 0000 1111 1111 1111
	}
	if (INDEX_DEPTH == 7) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 4096;		//binary number: 0000 0000 0001 0000 0000 0000
		BINARY_CODE[2] = 8192;		//binary number: 0000 0000 0010 0000 0000 0000
		BINARY_CODE[3] = 12288;		//binary number: 0000 0000 0011 0000 0000 0000
		BINARY_CODE[4] = 16383;		//binary number: 0000 0000 0011 1111 1111 1111
	}
	if (INDEX_DEPTH == 8) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 16384;		//binary number: 0000 0000 0100 0000 0000 0000
		BINARY_CODE[2] = 32768;		//binary number: 0000 0000 1000 0000 0000 0000
		BINARY_CODE[3] = 49152;		//binary number: 0000 0000 1100 0000 0000 0000
		BINARY_CODE[4] = 65535;		//binary number: 0000 0000 1111 1111 1111 1111
	}
	if (INDEX_DEPTH == 9) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 65536;		//binary number: 0000 0001 0000 0000 0000 0000
		BINARY_CODE[2] = 131072;	//binary number: 0000 0010 0000 0000 0000 0000
		BINARY_CODE[3] = 196608;	//binary number: 0000 0011 0000 0000 0000 0000
		BINARY_CODE[4] = 262143;	//binary number: 0000 0011 1111 1111 1111 1111
	}
	if (INDEX_DEPTH == 10) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 262144;	//binary number: 0000 0100 0000 0000 0000 0000
		BINARY_CODE[2] = 524288;	//binary number: 0000 1000 0000 0000 0000 0000
		BINARY_CODE[3] = 786432;	//binary number: 0000 1100 0000 0000 0000 0000
		BINARY_CODE[4] = 1048575;	//binary number: 0000 1111 1111 1111 1111 1111
	}
	if (INDEX_DEPTH == 11) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 1048576;	//binary number: 0001 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 2097152;	//binary number: 0010 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 3145728;	//binary number: 0011 0000 0000 0000 0000 0000
		BINARY_CODE[4] = 4194303;	//binary number: 0011 1111 1111 1111 1111 1111
	}
	if (INDEX_DEPTH == 12) {
		BINARY_CODE[0] = 0;			//binary number: 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 4194304;	//binary number: 0100 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 8388608;	//binary number: 1000 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 12582912;	//binary number: 1100 0000 0000 0000 0000 0000
		BINARY_CODE[4] = 16777215;	//binary number: 1111 1111 1111 1111 1111 1111
	}
	if (INDEX_DEPTH == 13) {
		BINARY_CODE[0] = 0;             //binary number: 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 16777216;      //binary number: 0001 0000 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 33554432;      //binary number: 0010 0000 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 50331648;      //binary number: 0011 0000 0000 0000 0000 0000 0000
		BINARY_CODE[4] = 67108863;		//binary number: 0011 1111 1111 1111 1111 1111 1111
	}
	if (INDEX_DEPTH == 14) {
		BINARY_CODE[0] = 0;             //binary number: 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 67108864;      //binary number: 0100 0000 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 134217728;     //binary number: 1000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 201326592;     //binary number: 1100 0000 0000 0000 0000 0000 0000
		BINARY_CODE[4] = 268435455;		//binary number: 1111 1111 1111 1111 1111 1111 1111
	}
	if (INDEX_DEPTH == 15) {
		BINARY_CODE[0] = 0;             //binary number: 0000 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[1] = 268435456;     //binary number: 0001 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[2] = 536870912;     //binary number: 0010 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[3] = 805306368;     //binary number: 0011 0000 0000 0000 0000 0000 0000 0000
		BINARY_CODE[4] = 268435456*4-1; //binary number: 0011 1111 1111 1111 1111 1111 1111 1111
	}

	return(0);
}

int init_genome_file() 
{
	char line[513], header = 0;

	if ((GENOME_FP = fopen(GENOME_FILE_NAME, "r")) == NULL) {
		fprintf(stderr, "ERROR : Couldn't open genome file %s\n", GENOME_FILE_NAME);
		exit(1);
	}

	do {
		if (fgets(line, 512, GENOME_FP) == NULL) {
			fprintf(stderr, "ERROR : Couldn't read input file %s\n", GENOME_FILE_NAME);
			exit(1);
		}
	} while (line[0] != '>');
	
	if (strchr(line, '\n') == NULL) header = 1;

	int chrdesclen = (strcspn(line, " \t\n") > CHR_DESC_LENGTH-1)? CHR_DESC_LENGTH-1: strcspn(line, " \t\n") - 1;
	if (chrdesclen > 0) strncpy(CHR_DESC, &line[1], chrdesclen);
		else {
			fprintf(stderr, "ERROR: A chromosome doesn't have a valid description!\n");
			exit(0);
		}
	CHR_DESC[chrdesclen] = '\0';
	strcpy(CHR_DESC_TMP, CHR_DESC);
	
	if (header == 1) {
		do {
			if (fgets(line, 512, GENOME_FP) == NULL) {
				fprintf(stderr, "ERROR : input file %s corrupt?\n", GENOME_FILE_NAME);
				exit(1);
			}
			if (strchr(line, '\n') != NULL) header = 0;
		} while (header != 0);
	}
	
	return 0;
}


int init_chr_index_file() 
{
	if ((CHR_INDEX_FP = fopen(CHR_INDEX_FILE_NAME, "w")) == NULL) {
		fprintf(stderr, "ERROR : Couldn't open index file %s\n", CHR_INDEX_FILE_NAME);
		exit(1);
	}

	return(0);
}

int init_mapindex_file()
{
        if ((MAPFWD_INDEX_FP = fopen(MAPFWD_INDEX_FILE_NAME, "w")) == NULL) {
                fprintf(stderr, "ERROR : Couldn't open mapindex file %s\n", MAPFWD_INDEX_FILE_NAME);
                exit(1);
        }

        return(0);
}


int init_meta_index_file() 
{
	if ((META_INDEX_FP = fopen(META_INDEX_FILE_NAME, "w")) == NULL) {
		fprintf(stderr, "ERROR : Couldn't open meta index file %s\n", META_INDEX_FILE_NAME);
		exit(1);
	}

	return(0);
}

int init_genome_out_file() 
{
	if ((GENOME_OUT_FP = fopen(GENOME_OUT_FILE_NAME, "w")) == NULL) {
		fprintf(stderr, "ERROR : Couldn't open meta index file %s\n", GENOME_OUT_FILE_NAME);
		exit(1);
	}

	return(0);
}

/*int init_stat_file()
{
	if ((OCC_FP = fopen(OCC_FILE_NAME, "w")) == NULL) {
		fprintf(stderr, "ERROR : Couldn't open stat file %s\n", OCC_FILE_NAME);
		exit(1);
	}

	return(0);
}*/

int init_mem_master() 
{
	if ( (MEM_MGR = (STORAGE *) malloc(sizeof(STORAGE))) == NULL) {
		fprintf(stderr, "ERROR : not enough memory for memory manager\n");
		exit(1);
	}

	MEM_MGR->curr_num = 0;
	MEM_MGR->nuggets = NULL;

	return 0;
}

