/* written by Jonas Behr, Regina Bohnert and Gunnar Raetsch, FML Tuebingen, Germany, 2010 */

#ifndef __GET_READS_DIRECT_H__
#define __GET_READS_DIRECT_H__

#include <vector>
  using std::vector;
#include "read.h"

static int g_flag_on = 0, g_flag_off = 0;

static int subsample = 100; 
static int collapse = 0;

int get_reads_from_bam(char* filename, char* region, vector<CRead*>* reads, char strand, int lsubsample);

#endif
