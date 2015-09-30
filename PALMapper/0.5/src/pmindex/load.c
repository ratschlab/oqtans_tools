// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "pmindex.h"

int load_chr();
int alloc_chr_seq_buffer();
int load_chr_sequence();


int load_chromosomes()
{
	unsigned int eof = 0;
	unsigned int chr = 0;
	POS p;
	
	BLOCK = 0;


	while (!eof) {

		if (VERBOSE) { printf("Start chromosome nb. %d\n", chr+1); }

		eof = load_chr();

		p.chr = chr;
		p.pos = 0;
		BLOCK_TABLE[BLOCK] = p;
		POSITION = 0;

		if (VERBOSE) { printf("\tLength %d\n", CHR_LENGTH); }


		index_chromosome(chr);

		write_chr_desc(chr);

		//if (!eof) {
		//	dealloc_chr();
		//}
		
		++BLOCK;
		++chr;

	}

	//write_index();

	if (VERBOSE) { printf("Write meta data to file ..."); fflush(stdout); }

	write_meta_index(chr);

	return(0);
}


int load_chr()
{
	char eof = 0;

	alloc_chr_seq_buffer();
	eof = load_chr_sequence();

	return(eof);
}

int load_chr_sequence()
{
	char line[513];

	char eof = 0; // end of file
	char eos = 0; // end of sequence
	unsigned int i, linelen, chrdesclen;
	
	CHR_LENGTH = 0;

	if (VERBOSE) { printf("\tReading chromosome into index ..."); }

	while (eof == 0 && eos == 0) {
		if (fgets(line, 512, GENOME_FP) == NULL) {
			eof = 1;
			eos = 1;
		}
		else {
			linelen = strcspn(line, "> \n\t");
			if (linelen != 0) {
				for (i=0; i != linelen; i++) {
					if (line[i]=='A' || line[i]=='a' || line[i]=='C' || line[i]=='c' ||
					    line[i]=='G' || line[i]=='g' || line[i]=='T' || line[i]=='t' || 
					    line[i]=='N' || line[i]=='n' || line[i]=='R' || line[i]=='r' ||
					    line[i]=='Y' || line[i]=='y' || line[i]=='M' || line[i]=='m' || 
					    line[i]=='K' || line[i]=='k' || line[i]=='W' || line[i]=='w' || 
					    line[i]=='S' || line[i]=='s' || line[i]=='B' || line[i]=='b' ||
					    line[i]=='D' || line[i]=='d' || line[i]=='H' || line[i]=='h' ||
					    line[i]=='V' || line[i]=='v') {
						CHR_SEQ[CHR_LENGTH] = toupper(line[i]);
						CHR_LENGTH++;
					}
					else {
						fprintf(stderr,"ERROR: Character '%c' encountered in chromosome '%s'! Only IUPAC-code is accepted!\n", line[i], CHR_DESC_TMP);
						exit(0);
					} 
				}
			}
			else {
				if (line[0] == '>') eos = 1;
			}
		}
	}
	
	CHR_SEQ[CHR_LENGTH] = '\0';

	if (CHR_LENGTH > LONGEST_CHROMOSOME) {
		LONGEST_CHROMOSOME = CHR_LENGTH;
	}

	strcpy(CHR_DESC, CHR_DESC_TMP);
	//if (strlen(CHR_DESC_TMP) != 0) CHR_DESC[strlen(CHR_DESC_TMP)] = '\0';

	if (eof == 0) {
		chrdesclen = (strcspn(line, " \t\n") > CHR_DESC_LENGTH-1)? CHR_DESC_LENGTH-1: strcspn(line, " \t\n") - 1;
		if (chrdesclen > 0) strncpy(CHR_DESC_TMP, &line[1], chrdesclen);
		else {
			fprintf(stderr, "ERROR: A chromosome doesn't have a valid description!\n");
			exit(0);
		}
		CHR_DESC_TMP[chrdesclen] = '\0';

		// to the next sequence (description can be longer than 512 chars)
		while (strchr(line,'\n') == NULL) {
			if (fgets(line, 512, GENOME_FP) == NULL) {
				fprintf(stderr, "ERROR: Missing sequence of chromosome/contig '%s'!\n", CHR_DESC);
				exit(0);
			}
		}
	}
	
	//CHR_DESC[strcspn(CHR_DESC, " \t")] = '\0';

	if (VERBOSE) { printf("... done\n"); printf("\tChromosome description: %s\n", CHR_DESC); }

	return(eof);
}


int alloc_chr_seq_buffer() 
{
	char line[513], header = 0;
	unsigned int fp, linelen;
	unsigned int l = 0;
	
	free(CHR_SEQ);

	fp = ftell(GENOME_FP);
	while (fgets(line, 512, GENOME_FP) != NULL) {
		
		if (header == 1 && strchr(line, '\n') != NULL) break;
		
		if (l == 0 && line[0] == '>') {
			fprintf(stderr, "ERROR: Missing sequence in reference input file!\n");
			exit(0);
		}
		
		linelen = strcspn(line, "> \n\t");
		
		if (header != 1) {
			if (linelen != 0) {
				if (line[linelen] == ' ' || line[linelen] == '\t') {
					fprintf(stderr, "ERROR: white space character unequal to newline found in genome input file '%s' in chromosome '%s'!\n", GENOME_FILE_NAME, CHR_DESC_TMP);
					exit(0);
				}
				l += linelen;
			}
			else {
				if (line[0] == '>') {
					if (strchr(line, '\n') == NULL) header = 1;
					if (strcspn(line+1, " \t\n") == 0) {
						fprintf(stderr, "ERROR: A chromosome doesn't have a description!\n");
						exit(0);
					}
					if (header == 0) break;
				}
			}
		}
	}

	if ((CHR_SEQ = (char *) malloc ((l + 1) * sizeof(char))) == NULL) {
		fprintf(stderr, "ERROR : couldn't allocate memory for index sequence\n");
		exit(1);
	}

	if((fseek(GENOME_FP, fp, SEEK_SET)) != 0) {
		fprintf(stderr, "ERROR: unable to move file pointer\n");
		exit(1);
	}

	return 0;
}
