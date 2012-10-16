#include <mex.h>
#include <stdlib.h>
char *get_string(const mxArray *prhs)
{
	char *buf;
	int buflen;

	if (!prhs)
		mexErrMsgTxt("get_string called with NULL pointer arg");
	if (!mxIsChar(prhs))
		mexErrMsgTxt("input is not a string");
	if (mxGetM(prhs) != 1)
		mexErrMsgTxt("input is not a row vector");
	buflen = mxGetN(prhs) + 1;
	buf = malloc(buflen);
	/* copy the string from prhs into buf and add terminating NULL char */
	if (mxGetString(prhs, buf, buflen))
		mexErrMsgTxt("not enough space");
	return buf;
}
