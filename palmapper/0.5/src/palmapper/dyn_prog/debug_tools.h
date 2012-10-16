#ifndef __DEBUG_TOOLS__
#define __DEBUG_TOOLS__

#include <cstdio>

void fassert(bool exp,int line, char* file);

#endif // __DEBUG_TOOLS__

#define FA(expr) (fassert(expr,__LINE__,__FILE__))

