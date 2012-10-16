// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "pmindex.h"

int usage() {
#ifdef GM
	printf("\ngmindex v%s\n", VERSION);
#else
	printf("\npmindex v%s\n", VERSION);
#endif
	printf("developed by Korbinian Schneeberger, Stephan Ossowski and Joerg Hagmann\n");
	printf("Max Planck Institute for Developmental Biology, Tübingen, Germany, 2008\n\n");
#ifdef GM
	printf("USAGE: gmindex [options]\n");
#else
        printf("USAGE: pmindex [options]\n");
#endif
	printf("\n");
	printf("mandatory:\n");
	printf(" -i STRING  input fastafile\n");
	//printf(" -x STRING  index filename\n");
	//printf(" -m STRING  map forward index filename\n");	
	//printf(" -c STRING  map reverse index filename\n");	
	//printf(" -t STRING  meta index filename\n");
	printf("\n");
	printf("optional:\n");
	printf(" -s INT     seed length, range 5 to 13 (default: 12)\n");
	printf(" -v         verbose (default: silent)\n");
	printf("\n");

	return 0;
}
