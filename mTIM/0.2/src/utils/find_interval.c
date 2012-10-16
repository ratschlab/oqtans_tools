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
#include <assert.h>
#include <mex.h>

#define MAXLINE 255
/*#define DEBUG 1*/

extern char *get_string(const mxArray *prhs);
extern int mmap_file(const char *filename, int open_mode, void **map, off_t *size) ;

/* returns 0 if an empty interval was found else the interval size, 
   lindex and rindex are only set if return value is positive */ 
int find_interval(unsigned *pos_map, off_t num_entries,
		   unsigned begin, unsigned end, unsigned *lindex, unsigned *rindex)
{
	unsigned val, mid, left = 0, right = num_entries - 1;
	assert(begin <= end);
	assert(num_entries > 0);
	if (pos_map[right] < begin || pos_map[left] > end) {
#ifdef DEBUG	
	  printf("ret 0\n\n");
#endif
	  rindex = 0 ;
	  lindex = 0 ;
          return 0; /* pos_map cannot possibly contain any values between begin and end */
	}

	for (;;) {
	  mid = (right + left) / 2; /*round down */
	  val = pos_map[mid];
#ifdef DEBUG	
	  printf("left %u, right %u, mid %u, val %u\n", left, right, mid, val);
#endif
	  if (val >= begin) {
	    if (right == mid)
	      break;
	    right = mid;
	  } else {
	    if (left == mid)
	        break;
	    left = mid;
	  }
	}
	*lindex = left;
	if (pos_map[left] < begin)
	  (*lindex)++;
#ifdef DEBUG	
	printf("left %u, right %u, mid %u, val %u\n", left, right, mid, val);
	printf("pos_map[left] %u, pos_map[left+1] %u, pos_map[left+2] %u\n",
	       pos_map[left], pos_map[left+1], pos_map[left+2]);
	printf("lindex %u, pos_map[lindex] %u\n", *lindex, pos_map[*lindex]);
#endif
	left = 0;
	right = num_entries - 1;
	for (;;) {
  	    mid = (right + left + 1) / 2; /*round up */
	    val = pos_map[mid];
#ifdef DEBUG	
	  printf("left %u, right %u, mid %u, val %u\n", left, right, mid, val);
#endif
	    if (val <= end) {
	        if (left == mid)
		    break;
		left = mid;
	    } else {
	      if (right == mid)
		  break;
	      right = mid;
	    }
	}
	*rindex = mid;
	if (pos_map[right] > end)
	  (*rindex)--;

#ifdef DEBUG	
	printf("left %u, right %u, mid %u, val %u\n", left, right, mid, val);
	printf("pos_map[right] %u, pos_map[right-1] %u, pos_map[right-2] %u\n", 
	       pos_map[right], pos_map[right-1], pos_map[right-2]);
	printf("rindex %u, pos_map[rindex] %u\n", *rindex, pos_map[*rindex]);
	printf("ret %i\n\n", *rindex - *lindex + 1);
#endif
	if (!(*rindex+1 >= *lindex))
	  {
	    fprintf(stderr, "assertion failed: *rindex + 1 = %i + 1 !>=! *lindex=%i\n", *rindex, *lindex) ;
	    mexErrMsgTxt("Assertion failed") ;
	  }
	/*if (!(*rindex >= *lindex))
	 *rindex=*lindex ;*/
        return *rindex - *lindex + 1;
}
