// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "pmindex.h"

// Deprecated:
int dealloc_chr() 
{

//	################ Delete information about positions in bins ######################
	int i, j;
	NUGGET *nugget, *next_nugget;
	
	for (j=0; j<NUM_USED_SLOTS; j++) {
                i = USED_SLOTS[j];

		if(INDEX[i] != NULL) {
			INDEX[i]->num_pos = 0;
			INDEX[i]->bin_ext = 0;
			INDEX[i]->last_bin_ext = 0;
		}
	}

	NUM_USED_SLOTS = 0;

//	################ Dealloc the extended bins ######################
	nugget = MEM_MGR->nuggets;
	while (nugget != NULL) {
		next_nugget = nugget->next;
		free(nugget);
		nugget = next_nugget;
	}

	MEM_MGR->curr_num = 0;
	MEM_MGR->nuggets = NULL;

	return 0;
}

int alloc_bin(int slot) {

	if (((INDEX[slot] = (BIN *) malloc (sizeof(BIN))) == NULL)) {
		fprintf(stderr, "ERROR : couldn't allocate memory for a storage bin\n");
		exit(1);
	}
	
	INDEX[slot]->num_pos = 0;
	INDEX[slot]->bin_ext = 0;
	INDEX[slot]->last_bin_ext = 0;

	return 0;
}

BIN_EXT *alloc_bin_ext() 
{ 

	NUGGET *tmp;
	BIN_EXT *bin_ext;
	
	if (MEM_MGR->nuggets == NULL || MEM_MGR->curr_num == (BIN_EXT_PER_NUGGET-1)) { 
		tmp = MEM_MGR->nuggets; 

		if ((MEM_MGR->nuggets = (NUGGET *) malloc (sizeof(NUGGET))) == 0) {
			fprintf (stderr, "ERROR : could not allocate memory for nugget storage\n");
			exit(1);
		}
		
		(MEM_MGR->nuggets)->next = tmp;
		MEM_MGR->curr_num = -1; 
	}

	MEM_MGR->curr_num++;

	bin_ext = &((MEM_MGR->nuggets)->buffer[MEM_MGR->curr_num]);
	bin_ext->bin_ext = 0;

	return(bin_ext);
}

void alloc_blocktable()
{
	BLOCK_TABLE = (POS *) malloc (BLOCK_TABLE_SIZE * sizeof(POS));
	
}
