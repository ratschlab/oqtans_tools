/*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3 of the License, or
*  (at your option) any later version.
*
*   Written (W) 2010-2011 Jonas Behr
*   Copyright (C) 2010-2011 Max Planck Society
*/


#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/mman.h>
#include <time.h>
#include <math.h>
#include <limits>
#include <mex.h>
#include <assert.h>
#include <vector>
  using std::vector;
#include <algorithm>
  using std::sort;
  using std::min;
  using std::max;

typedef struct {
    int start;
    int stop;
	int idx;
	int set_id;
} interval_t;

bool compare (interval_t i, interval_t j) 
{ 
	return (i.start<j.start); 
}

bool overlaps(interval_t a, interval_t b)
{
	int v = min(a.stop,b.stop) - max(a.start,b.start) + 1;
	return (v >= 1);
}
bool leftOf(interval_t a, interval_t b)
{
	return (a.stop < b.start);
}

void scan(interval_t f, vector<interval_t>* Wf, interval_t g, vector<interval_t>* Wg, vector<int>* overlap)
{
	vector<interval_t>::iterator i;
	i=Wg->begin();
	while (i<Wg->end())
	{
		interval_t g2 = *i;
		if (leftOf(g2,f))
		{
			Wg->erase(i);// inefficient if Wg is large
			// this moves all elements, therefore i is not incremented
		}
		else if (overlaps(g2,f))
		{
			if (g2.set_id==1)
			{
				//printf("overlap: [%i | %i, %i] [%i | %i, %i]\n", g2.idx, g2.start, g2.stop, f.idx, f.start, f.stop);
				overlap->push_back(g2.idx);
				overlap->push_back(f.idx);
			}
			else if (f.set_id==1)
			{
				//printf("overlap: [%i | %i, %i] [%i | %i, %i]\n", f.idx, f.start, f.stop, g2.idx, g2.start, g2.stop);
				overlap->push_back(f.idx);
				overlap->push_back(g2.idx);
			}	
			i++;
		}
		else
		{
			printf("never happens??\n");
			i++;
		}
	}
	if (!leftOf(f, g))
	{
		Wf->push_back(f);
		//printf("push: [%i, %i] size:%i\n", f.start, f.stop, Wf->size());
	}
}

/*
 * prhs[0] first interval set starts
 * prhs[1] first interval set stops
 * prhs[2] second interval set starts
 * prhs[3] second interval set stops
 *
 * return:
 * plhs[0] one based index in first interval set overlapping with a interval in the second set
 * plhs[1] corresponding index in the second set
 *
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs!=4)
		mexErrMsgTxt("Expected 4 arguments: starts1, stops1, starts2, stops2 \n");
	if (nlhs!=2)
		mexErrMsgTxt("Expected 2 output arguments \n");

	int num_intervals1 = mxGetNumberOfElements(prhs[0]);
	assert(num_intervals1 == mxGetNumberOfElements(prhs[1]));
	int num_intervals2 = mxGetNumberOfElements(prhs[2]);
	assert(num_intervals2 == mxGetNumberOfElements(prhs[3]));
	
	//printf("num_intervals1: %i\n", num_intervals1);
	//printf("num_intervals2: %i\n", num_intervals2);

	double* starts1 = mxGetPr(prhs[0]);
	double* stops1  = mxGetPr(prhs[1]);
	double* starts2 = mxGetPr(prhs[2]);
	double* stops2  = mxGetPr(prhs[3]);
	
	vector<interval_t> intervals1;
	for (int i=0; i<num_intervals1; i++)
	{
		interval_t interval;
		interval.start = starts1[i];
		interval.stop = stops1[i];
		interval.set_id = 1;
		interval.idx = i+1;
		intervals1.push_back(interval);
		//printf("int1: [%i, %i] \n",intervals1[i].start, intervals1[i].stop);
	}
	interval_t i;
	i.start = std::numeric_limits<int>::max();
	i.stop = std::numeric_limits<int>::max();
	i.set_id = std::numeric_limits<int>::max();
	i.idx = std::numeric_limits<int>::max();
	intervals1.push_back(i);

	//printf("num_intervals1: %i\n", intervals1.size());
	vector<interval_t> intervals2;
	for (int i=0; i<num_intervals2; i++)
	{
		interval_t interval;
		interval.start = starts2[i];
		interval.stop = stops2[i];
		interval.set_id = 2;
		interval.idx = i+1;
		intervals2.push_back(interval);
		//printf("int2: [%i, %i] \n",intervals2[i].start, intervals2[i].stop);
	}
	intervals2.push_back(i);
	//printf("num_intervals2: %i\n", intervals2.size());

	sort(intervals1.begin(), intervals1.end(), compare);	
	sort(intervals2.begin(), intervals2.end(), compare);	


	vector<int> overlap;
	vector<interval_t> Wx;
	vector<interval_t> Wy;
	vector<interval_t>::iterator x = intervals1.begin();
	vector<interval_t>::iterator y = intervals2.begin();
	while (x<intervals1.end() && y<intervals2.end())
	{
		//vector<interval_t>::iterator x;
		//vector<interval_t>::iterator y;
		//if (it1>intervals1.end())
		//	x = inf_interval();
		//else
		//	x = it1;
		//if (it2>intervals2.end())
		//	y = inf_interval();
		//else
		//	y=it2;

		if (x->start <= y->start)
		{
				scan(*x, &Wx, *y, &Wy, &overlap);
				x++;
		}
		else
		{
			if (y<=intervals2.end())
			{
				scan(*y, &Wy, *x, &Wx, &overlap);
				y++;
			}
		}
	}

	plhs[0] = mxCreateDoubleMatrix(1, overlap.size()/2, mxREAL);
    double* idx1 = mxGetPr(plhs[0]);

	plhs[1] = mxCreateDoubleMatrix(1, overlap.size()/2, mxREAL);
    double* idx2 = mxGetPr(plhs[1]);

	for (int i=0; i<overlap.size(); i+=2)
	{
		//printf("ov: %i [%i, %i] \n", i, overlap[i], overlap[i+1]);
		idx1[i/2] = (double) overlap[i];
		idx2[i/2] = (double) overlap[i+1];
	}
}





