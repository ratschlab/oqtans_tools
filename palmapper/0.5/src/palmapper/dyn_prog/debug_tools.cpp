// Authors: Uta Schulze, Cheng Soon Ong, Fabio De Bona, Gunnar Raetsch, Geraldine Jean
// Copyright (C) 2006-2010 by Friedrich Miescher Laboratory, Tuebingen, Germany

#include <cstdio>
#include <cstdlib>
#include "debug_tools.h"
using namespace std;

void fassert(bool exp,int line, char* file) {
   if (! exp) {
      printf("invalid fassert at line %d in file %s!\n",line,file);
      exit(EXIT_FAILURE);
   }
}
