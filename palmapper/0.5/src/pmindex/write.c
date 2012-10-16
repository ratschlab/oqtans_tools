// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "pmindex.h"

int write_chr_desc(unsigned int chr) 
{
	if (fwrite(&chr, sizeof(unsigned int), 1, CHR_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                exit(0);
        }
        if (fwrite(&CHR_LENGTH, sizeof(unsigned int), 1, CHR_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                exit(0);
        }
        if (fwrite(&CHR_DESC[0], sizeof(char), CHR_DESC_LENGTH, CHR_INDEX_FP) == 0) {
                fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                exit(0);
        }       //@TODO CHR_LENGTH instead of CHR_DESC_LENGTH

	

	return 0;
}

/**
int write_index()
{

        unsigned int j, i, num;
        int minus_i;
        BIN *bin;
        BIN_EXT **bin_ext;

	for (j=0; j<NUM_USED_SLOTS; j++) {
		i = USED_SLOTS[j];
		
		if (INDEX[i] != NULL && INDEX[i]->num_pos != 0) {

			bin = INDEX[i];
			num = bin->num_pos;

			//HEADER
			if (fwrite(&i, sizeof(int), 1, INDEX_FP) == 0) {	// slot
				fprintf(stderr, "ERROR: cant access harddisc for index file\n"); 
			}

			if (fwrite(&num, sizeof(int), 1, INDEX_FP) == 0) {	// nr of positions in this slot
				fprintf(stderr, "ERROR: cant access harddisc for index file\n");
			}

			//CONTENT ( = positions)
			if (num <= BIN_SIZE) {

				if (fwrite(&bin->ids[0], sizeof(ID), num, MAPFWD_INDEX_FP) == 0) {
                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                        exit(0);
                                }
			
			} else { 

				if (fwrite(&(bin->ids[0]), sizeof(ID), BIN_SIZE, MAPFWD_INDEX_FP) == 0) {
                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                        exit(0);
                                }
				num -= BIN_SIZE;

				bin_ext = &(bin->bin_ext);

				while (*bin_ext != NULL) {
					if (num > BIN_SIZE_EXT) {  

						if (fwrite(&((*bin_ext)->ids[0]), sizeof(ID), BIN_SIZE_EXT, MAPFWD_INDEX_FP) == 0) {
                                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                                        exit(0);
                                                }
						num -= BIN_SIZE_EXT;
						bin_ext = &((*bin_ext)->bin_ext);

					} else { //write out last entries and finish

						if (fwrite(&((*bin_ext)->ids[0]), sizeof(ID), num, MAPFWD_INDEX_FP) == 0) {
                                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                                        exit(0);
                                                }
						*bin_ext = NULL;

					}
				}
			}
		}

		
		// and the same for the reverse index
		if (BUILD_REVERSE_INDEX && INDEX_REV[i] != NULL && INDEX_REV[i]->num_pos != 0) {
					
			bin = INDEX_REV[i];
			num = bin->num_pos;

			//HEADER
			//slot
			if (i == 0) minus_i = -2147483647;
				else minus_i = -i;
			
			// slot writting	
			if (fwrite(&minus_i, sizeof(int), 1, INDEX_FP) == 0) {
				fprintf(stderr, "ERROR: cant access harddisc for index file\n"); 
				exit(0);
			}

			// number of entries
			if (fwrite(&num, sizeof(int), 1, INDEX_FP) == 0) {
				fprintf(stderr, "ERROR: cant access harddisc for index file\n");
				exit(0);
			}

			//CONTENT
			if (num <= BIN_SIZE) {

				if (fwrite(&bin->ids[0], sizeof(ID), num, MAPREV_INDEX_FP)  == 0) {
                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                        exit(0);
                                }

			} else { 

				if (fwrite(&(bin->ids[0]), sizeof(ID), BIN_SIZE, MAPREV_INDEX_FP) == 0) {
                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                        exit(0);
                                }
				num -= BIN_SIZE;

				bin_ext = &(bin->bin_ext);

				while (*bin_ext != NULL) {
					if (num > BIN_SIZE_EXT) {  
						if (fwrite(&((*bin_ext)->ids[0]), sizeof(ID), BIN_SIZE_EXT, MAPREV_INDEX_FP) == 0) {
                                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                                        exit(0);
                                                }
						num -= BIN_SIZE_EXT;
						bin_ext = &((*bin_ext)->bin_ext);

					} else { //write out last entries and finish
						if (fwrite(&((*bin_ext)->ids[0]), sizeof(ID), num, MAPREV_INDEX_FP) == 0) {
                                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                                        exit(0);
                                                }
						*bin_ext = NULL;					

					}
				}
			}
		}
	}
	
	if (VERBOSE) printf("... done\n");

	return(0);
}
*/

int write_meta_index(unsigned int num_chr) 
{
	BIN *bin = 0;
        BIN_EXT **bin_ext;

	unsigned int i, num;


	// write meta information
	char dummy=0;	//downcompatibility to previous GM versions:
	if (fwrite(&dummy, sizeof(char), 1, META_INDEX_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisc for meta index file\n");
	}
	if (fwrite(&INDEX_DEPTH, sizeof(int), 1, META_INDEX_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisc for meta index file\n"); 
	}
	if (fwrite(&num_chr, sizeof(unsigned int), 1, META_INDEX_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisc for meta index file\n"); 
	}
	if (fwrite(&POSITION_COUNTER, sizeof(unsigned int), 1, META_INDEX_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisc for meta index file\n"); 
	}
	if (fwrite(&LONGEST_CHROMOSOME, sizeof(unsigned int), 1, META_INDEX_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisc for meta index file\n"); 
	}


	// write block table
	if (fwrite(&BLOCK, sizeof(unsigned int), 1, META_INDEX_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisc for meta index file\n"); 
	}	
	if (fwrite(BLOCK_TABLE, sizeof(POS), BLOCK, META_INDEX_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisc for meta index file\n"); 
	}
	

	
	// write bins: 
	unsigned long int maxnr = 0;

 	for (i=0; i<INDEX_SIZE; i++) {

		if ( INDEX[i] != NULL ) {	

			bin = INDEX[i];
                        num = bin->num_pos;

			//Slot
			if (fwrite(&i, sizeof(int), 1, META_INDEX_FP) == 0) {
				fprintf(stderr, "ERROR: cant access harddisc for meta index file\n"); 
			}

			//Number
			if (fwrite(&num, sizeof(int), 1, META_INDEX_FP) == 0) {
				fprintf(stderr, "ERROR: cant access harddisc for index file\n");
			}
			maxnr = (num > maxnr)? num : maxnr;

			//CONTENT
                        if (num <= BIN_SIZE) {

                                if (fwrite(&bin->ids[0], sizeof(ID), num, MAPFWD_INDEX_FP) == 0) {
                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                        exit(0);
                                }

                        } else {

                                if (fwrite(&(bin->ids[0]), sizeof(ID), BIN_SIZE, MAPFWD_INDEX_FP) == 0) {
                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                        exit(0);
                                }
                                num -= BIN_SIZE;

                                bin_ext = &(bin->bin_ext);

                                while (*bin_ext != NULL) {
                                        if (num > BIN_SIZE_EXT) {

                                                if (fwrite(&((*bin_ext)->ids[0]), sizeof(ID), BIN_SIZE_EXT, MAPFWD_INDEX_FP) == 0) {
                                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                                        exit(0);
                                                }
                                                num -= BIN_SIZE_EXT;
                                                bin_ext = &((*bin_ext)->bin_ext);

                                        } else { //write out last entries and finish

                                                if (fwrite(&((*bin_ext)->ids[0]), sizeof(ID), num, MAPFWD_INDEX_FP) == 0) {
                                                        fprintf(stderr, "ERROR: cant access harddisc for index file\n");
                                                        exit(0);
                                                }
                                                *bin_ext = NULL;

                                        }
                                }
                        }
		}

	}

	if (VERBOSE) printf("... done\n");
	
	return(0);
}

/**
int write_chr(unsigned int chr)
{
	if (VERBOSE) printf("writing chromosome %d\n",chr+1);
	if (fwrite(CHR_SEQ, sizeof(char), CHR_LENGTH, GENOME_OUT_FP) == 0) {
		fprintf(stderr, "ERROR: cant access harddisk for genome output file\n");
	}
	
	fclose(GENOME_OUT_FP);
	
	return 0;	
}
*/

