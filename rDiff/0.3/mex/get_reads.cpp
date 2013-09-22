/* written by Jonas Behr, Regina Bohnert, Philipp Drewe and Gunnar Raetsch, FML Tuebingen, Germany, 2012 */

#include <stdio.h>
#include <signal.h>
#include <mex.h>
#include <vector>
	using std::vector;
#include "get_reads_direct.h"
#include "mex_input.h"
#include "read.h"
#include <sys/types.h>
#define MAXLINE 10000

/*
 * input: 
 * 1 bam file
 * 2 chromosome
 * 3 region start (1-based index)
 * 4 region end (1-based index)
 * 5 strand (either '+' or '-' or '0')
 * [6] collapse flag: if true the reads are collapsed to a coverage track
 * [7] subsample percentage: percentage of reads to be subsampled (in per mill)
 * [8] intron length filter
 * [9] exon length filter
 * [10] mismatch filter
 *
 * output: 
 * 1 coverage
 * [2] intron cell array
 *
 * example call: 
 * [cov introns] = get_reads('polyA_left_I+_el15_mm1_spliced.bam', 'I', 10000, 12000, '-', 1, 30);
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  if (nrhs<5 || nrhs>10 || (nlhs!=1 && nlhs!=2)) {
    fprintf(stderr, "usage: [x [introns]] = get_reads(fname, chr, start, end, strand, [collapse], [subsample], [max intron length], [min exonlength], [max mismatches]);\n");
    return; 
  }
  
  /* obligatory arguments
   * **********************/
  char *fname = get_string(prhs[0]);
	//fprintf(stdout, "arg1: %s\n", fname);
  char *chr = get_string(prhs[1]);
	//fprintf(stdout, "arg2: %s\n", chr);
  int from_pos = get_int(prhs[2]);
	//fprintf(stdout, "arg3: %d\n", from_pos);
  int to_pos = get_int(prhs[3]);
	//fprintf(stdout, "arg4: %d\n", to_pos);
  char *strand = get_string(prhs[4]);
	//fprintf(stdout, "arg5: %s\n", strand);

  if (from_pos>to_pos)
     mexErrMsgTxt("Start (arg 3) must be <= end (arg 4)\n");

  if (strand[0]!='+' && strand[0]!='-' && strand[0]!='0') 
    mexErrMsgTxt("Unknown strand (arg 5): either + or - or 0");

  /* optional arguments
   * ******************/	
  int collapse = 0;
  if (nrhs>=6)
    collapse = get_int(prhs[5]);
  
  int subsample = 1000;	
  if (nrhs>=7)
    subsample = get_int(prhs[6]);
	  
  int intron_len_filter = 1e9;
  if (nrhs>=8)
    intron_len_filter = get_int(prhs[7]);

  int exon_len_filter = -1;
  if (nrhs>=9)
    exon_len_filter = get_int(prhs[8]);

  int filter_mismatch = 1e9;
  if (nrhs>=10)
    filter_mismatch = get_int(prhs[9]);

  /* call function to get reads
   * **************************/	
  char region[MAXLINE];
  sprintf(region, "%s:%i-%i", chr, from_pos, to_pos);

  vector<CRead*> all_reads;
  
  get_reads_from_bam(fname, region, &all_reads, strand[0], subsample);
 
	/* filter reads
	* **************/	
	vector<CRead*> reads;
	for (int i=0; i<all_reads.size(); i++) 
	{
		if (all_reads[i]->max_intron_len()<intron_len_filter && all_reads[i]->min_exon_len()>exon_len_filter && all_reads[i]->get_mismatches()<=filter_mismatch && all_reads[i]->multiple_alignment_index==0)
			reads.push_back(all_reads[i]);
	}

 
  /* prepare output
   * **************/	
  int num_rows = reads.size();
  int num_pos = to_pos-from_pos+1;
  
  // read coverages collapsed 
  if (collapse) {
    plhs[0] = mxCreateNumericMatrix(1, num_pos, mxUINT32_CLASS, mxREAL);
    uint32_t *mask_ret = (uint32_t*) mxGetData(plhs[0]);
    if (num_pos>0 && mask_ret==NULL)
      mexErrMsgTxt("Error allocating memory\n");
    for (int i=0; i<reads.size(); i++) 
	{
        reads[i]->get_coverage(from_pos, to_pos, mask_ret);
    }
  }
  // reads not collapsed
  else {
    uint32_t nzmax = 0; // maximal number of nonzero elements 
    int len = to_pos-from_pos+1;
    for (uint i=0; i<reads.size(); i++) {
	for (uint n = 0; n < reads[i]->block_starts.size(); n++) {
	  uint32_t from, to;
	  if (reads[i]->block_starts[n]+reads[i]->start_pos-from_pos >= 0)
	    from = reads[i]->block_starts[n]+reads[i]->start_pos-from_pos;
	  else
	    from = 0;
	  if (reads[i]->block_starts[n]+reads[i]->start_pos-from_pos+reads[i]->block_lengths[n] >= 0)
	    to = reads[i]->block_starts[n]+reads[i]->start_pos-from_pos+reads[i]->block_lengths[n];
	  else
	    to = 0;
	  for (int bp=from; bp<to&bp<len; bp++) {
	    nzmax++;
	  }
	}
    }
    // 1st row: row indices
    // 2nd row: column indices
    plhs[0] = mxCreateDoubleMatrix(2, nzmax, mxREAL);
    double *mask_ret = (double*) mxGetData(plhs[0]);
    if (nzmax>0 && mask_ret==NULL)
      mexErrMsgTxt("Error allocating memory\n");
    uint32_t mask_ret_c = 0; // counter
    for (uint i=0; i<reads.size(); i++) 
	{
		reads[i]->get_reads_sparse(from_pos, to_pos, mask_ret, mask_ret_c, i);
    }
    if (mask_ret_c!=2*nzmax)
      mexErrMsgTxt("Error filling index arrays for sparse matrix\n");
    
  }
  // introns
  if (nlhs==2) 
  {
    vector<int> intron_list;
	for (int i=0; i<reads.size(); i++) 
	{
		reads[i]->get_introns(&intron_list);
	}

    plhs[1] = mxCreateNumericMatrix(2, intron_list.size()/2, mxUINT32_CLASS, mxREAL);
	uint32_t *p_intron_list = (uint32_t*) mxGetData(plhs[1]);
	for (int p = 0; p<intron_list.size(); p++) 
	{	
		p_intron_list[p] = intron_list[p];
	}
	intron_list.clear();	
  }
  //// introns
  //if (nlhs==2) 
  //{
  //  // use x = [introns{:}]' in matlab to get the array of introns
  //  const int cnum_rows = num_rows;
  //  plhs[1] = mxCreateCellArray(1, &cnum_rows);
  //  mxArray *intron_cell = plhs[1];
  //  if (num_rows>0 && intron_cell==NULL)
  //    mexErrMsgTxt("Error allocating memory\n");
  //  vector<int> intron_list;
  //  for (int i=0; i<reads.size(); i++) 
  //  {
  //    reads[i]->get_introns(&intron_list);
  //    mxArray * mx_intron_list = mxCreateNumericMatrix(2, intron_list.size()/2, mxUINT32_CLASS, mxREAL);		
  //    uint32_t *p_intron_list = (uint32_t*) mxGetData(mx_intron_list);
  //    
  //    for (int p = 0; p<intron_list.size(); p++) 
  //    {
  //      p_intron_list[p] = intron_list[p];
  //    }
  //    mxSetCell(intron_cell, i, mx_intron_list);
  //    intron_list.clear();	
  //  }
  //}
	for (int i=0; i<all_reads.size(); i++)
		delete all_reads[i];

}

