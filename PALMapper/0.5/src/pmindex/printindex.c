// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "pmindex.h"

void printindex()
{
	int i;
	BIN* bin;

	printf("------%c--%c----\n",CHR_SEQ[CHR_LENGTH-2], CHR_SEQ[CHR_LENGTH-1]);
	for (i=0; i!=INDEX_SIZE; ++i) {
		if (*(INDEX+i) != NULL) {
			bin = *(INDEX+i);
			//printf("%i: %i %i\n", i, bin->positions[0].chr, bin->positions[0].pos);
		}
	}
	printf("------------\n");

}
