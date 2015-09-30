/*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3 of the License, or
*  (at your option) any later version.
*
*   Written (W) 2009-2011 Regina Bohnert
*   Copyright (C) 2009-2011 Max Planck Society
*/


#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <ctype.h>
#include <assert.h>

#include <vector>
  using std::vector;
#include <string>
  using std::string;
#include <algorithm>
  using std::find; 
  using std::min;

#include <mex.h>


char *get_string(const mxArray *prhs);

typedef unsigned int uint32_t;
typedef unsigned char uint8_t;

/*
 * [read_len num_reads] = get_bam_properties(fname, path_samtools, contig_name)
 *
 * -- input --
 * prhs[0] file name of paired reads in BAM format (sorted by read id)
 * prhs[1] path to samtools
 * prhs[2] contig name
 *
 * -- output --
 * plhs[0] length of read
 * plhs[1] number of unique reads
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // checks for the right number of arguments
  if (nrhs !=3 || nlhs > 2) {
    mexErrMsgTxt("number of input and output args should be 3 and 2\nUSAGE:\n     [read_len, num_reads] = get_bam_properties(fname, path_samtools, contig_name)\n");
    return;
  }

  signal(SIGCHLD, SIG_IGN); // avoid zombies

  // read input arguments
  char *fname = get_string(prhs[0]);
  char *path_samtools = get_string(prhs[1]);
  char *contig_name = get_string(prhs[2]);
  char command[10000];
  
  sprintf(command, "%s./samtools view %s %s 2>/dev/null", path_samtools, fname, contig_name);
  //printf("%s\n", command);

  // get number of unique reads
  int status;
  uint32_t num_unique_reads = 0;
  char command2[10000];
  sprintf(command2, "%s | cut -f 1 | sort -u | wc -l", command);
  FILE* fp = popen(command2, "r");
  if (fp == NULL) {
    mexErrMsgTxt("Error using popen\n");
  }
  int num_scans = 1;
  num_scans = fscanf(fp, "%d", &num_unique_reads);
  if (num_scans != 1) {
    rewind(fp);
    char ret[1000];
    fgets(ret, 1000, fp);
    fprintf(stdout, "%s", ret);
    mexErrMsgTxt("Could not determine number of reads\n");
  }
  status = pclose(fp);
  //printf("%i", num_unique_reads);
  
  // select reads for given positions and strand
  int num_rows_selected = min((int) num_unique_reads, 100);
  sprintf(command, "%s | head -n %i | cut -f 1-11", command, num_rows_selected);
  fp = popen(command, "r");
  if (fp == NULL) {
    mexErrMsgTxt("Error using popen\n");
  }
  /* SAM format
     1: read id, 2: flag, 3: reference name, 4: start (1-based, incl.), 5: mapping quality,
     6: CIGAR, 7: mate reference name, 8: mate start (1-based, incl.), 9: insert size, 10: read, 11: quality
     12+: additional tags  
  */
  uint32_t read_idx = 0, row_idx = 0, num_col = 0;
  uint32_t flag = 0, start_pos = 0, map_score = 0, mate_end_pos = 0, num_matches = 0, num_del = 0, num_ins = 0, ins_size = 0;
  char ri [1000], read_contig_name [1000], cg [1000], mate_read_id [1000], read [1000], read_qual [1000];
  string last_read_id;
  vector<uint32_t> block_lengths, block_starts;
  vector<string> read_ids;
  vector<string>::iterator it;
  
  uint32_t read_len = 0;
  bool empty_line = true;
  int num_rows = 0;
  while(empty_line && num_rows < num_rows_selected) {
    num_col = fscanf(fp, "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s", &ri, &flag, &read_contig_name, &start_pos, &map_score, &cg, &mate_read_id, &mate_end_pos, &ins_size, &read, &read_qual);
    if (num_col != 11) {
      mexErrMsgTxt("error reading SAM line\n");
    }
    
    string cigar = (string) cg;
    // ignore lines with reads w/o mapping information 
    if (start_pos == 0 || cigar.compare("*")==0) {
      continue;
    }
    // parse CIGAR
    uint last_c = 0;
    string last_str;
    num_matches = 0;
    char *end = NULL;
    uint32_t tmp_nm = 0, tmp_nd = 0, tmp_ni = 0;
    uint32_t last_block_start = 0, last_block_length = 0, last_intron_len = 0;
    block_lengths.clear(); block_starts.clear();
    
    for (uint c = 0; c < cigar.size(); c++) {
      switch (cigar[c]) {
      case 'M':
	last_str = cigar.substr(last_c, c-last_c);
	tmp_nm = strtoul(last_str.c_str(), &end, 10);
	if (*end != '\0')
	  mexErrMsgTxt("error: number of mismatches\n");
	end = NULL;
	last_block_length += tmp_nm;
	num_matches += tmp_nm;
	last_c = c + 1;
	break;
      case 'I':
	last_str = cigar.substr(last_c, c-last_c);
	tmp_ni = strtoul(last_str.c_str(), &end, 10);
	if (*end != '\0')
	  mexErrMsgTxt("error: number of insertions\n");
	end = NULL;
	num_ins += tmp_ni;
	last_c = c + 1;
	break;
      case 'D':
	last_str = cigar.substr(last_c, c-last_c);
	tmp_nd = strtoul(last_str.c_str(), &end, 10);
	if (*end != '\0')
	  mexErrMsgTxt("error: number of deletions\n");
	end = NULL;
	num_del += tmp_nd;
	last_block_length += tmp_nd;
	last_c = c + 1;
	break;
      case 'N':
	last_str = cigar.substr(last_c, c-last_c);
	last_intron_len = strtoul(last_str.c_str(), &end, 10);
	end = NULL;
	last_c = c + 1;
	break;
      case 'S':
	break;
      case 'H':
	break;
      case 'P':
	break;
      default:
	break;
      }
      if (cigar[c] == 'N' || c==cigar.size()-1) {
	block_starts.push_back(last_block_start);
	last_block_start = last_block_start + last_block_length + last_intron_len;
	last_intron_len = 0;
	block_lengths.push_back(last_block_length);
	last_block_length = 0;
      }
    }
    read_len = 0;
    for (uint n = 0; n < block_lengths.size(); n++) {
      read_len += block_lengths[n];
    }
    empty_line = false;
  } // end of stream parsing
	
  status = pclose(fp);
  
  if (empty_line) 
    mexErrMsgTxt("Could not determine read length\n");
  
  plhs[0] = mxCreateDoubleScalar((double) read_len);
  plhs[1] = mxCreateDoubleScalar((double) num_unique_reads);
  
  return;
}


char *get_string(const mxArray *prhs) {
  char *buf;
  int buflen;
  if (!prhs)
    mexErrMsgTxt("get_string called with NULL pointer arg");
  if (!mxIsChar(prhs))
    mexErrMsgTxt("input is not a string");
  if (mxGetM(prhs) != 1)
    mexErrMsgTxt("input is not a row vector");
  buflen = mxGetN(prhs) + 1;
  buf = (char*) malloc(buflen);
  /* copy the string from prhs into buf and add terminating NULL char */
  if (mxGetString(prhs, buf, buflen))
    mexErrMsgTxt("not enough space");
  return buf;
}
