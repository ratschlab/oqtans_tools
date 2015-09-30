/*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3 of the License, or
*  (at your option) any later version.
*
*   Written (W) 2010-2011 Jonas Behr, Regina Bohnert, Gunnar Raetsch
*   Copyright (C) 2010-2011 Max Planck Society
*/


#include <stdio.h>
#include <mex.h>
#include "mex_input.h"

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

bool get_bool(const mxArray *prhs)
{
	const int M = mxGetM(prhs);
	const int N = mxGetN(prhs);
	double *f = (double*) mxGetPr(prhs);

	if (!prhs)
		mexErrMsgTxt("Arg is NULL pointer");
	if (M != 1 || N != 1)
		mexErrMsgTxt("Arg is not a scalar");
	if (f[0] != 0)
		return true;
	return false;
}

int get_int(const mxArray *prhs)
{
	const int M = mxGetM(prhs);
	const int N = mxGetN(prhs);
	double *f = (double*) mxGetPr(prhs);

	if (!prhs)
		mexErrMsgTxt("Arg is NULL pointer");
	if (M != 1 || N != 1)
		mexErrMsgTxt("Arg is not a scalar");

	return (int) f[0];
}
