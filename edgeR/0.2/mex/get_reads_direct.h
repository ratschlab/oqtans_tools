/*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3 of the License, or
*  (at your option) any later version.
*
*   Written (W) 2010-2011 Jonas Behr, Regina Bohnert, Gunnar Raetsch
*   Copyright (C) 2010-2011 Max Planck Society
*/


#ifndef __GET_READS_DIRECT_H__
#define __GET_READS_DIRECT_H__

#include <vector>
  using std::vector;
#include "read.h"

//static int g_flag_on = 0, g_flag_off = 0;
static int left_flag_mask = strtol((char*) "0x40", 0, 0);
static int right_flag_mask = strtol((char*) "0x80", 0, 0);
static int reverse_flag_mask = strtol((char*) "0x10", 0, 0);

static int subsample = 1000; 
//static int collapse = 0;

int get_reads_from_bam(char* filename, char* region, vector<CRead*>* reads, char strand, int lsubsample);

#endif
