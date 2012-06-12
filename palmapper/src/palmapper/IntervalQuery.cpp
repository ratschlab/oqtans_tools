/*

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Written (W) 2007-2009 Fabio De Bona, Andre Noll, Gunnar Raetsch
Copyright (C) 2007-2009 Max-Planck-Society

*/

#include <cstdlib>
#include <cstdio>
using namespace std;

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
#include <pthread.h>

#include "IntervalQuery.h"

// this switch defines whether the float score maps should be converted to unsigned short int
// to save memory
#define CONVERT_SCORE_MAPS

#define MAXLINE 255
/* #define DEBUG 1 */

/*enum {E_DUMMY, E_SYNTAX};*/
enum {E_SYNTAX = 1, E_BAD_INTERVAL};

std::vector<std::string> IntervalQuery::mmap_fname_list ;
std::vector<void*> IntervalQuery::mmap_ptr_list ;
std::vector<int> IntervalQuery::mmap_size_list ;
clock_t IntervalQuery::total_time ;
pthread_mutex_t IntervalQuery::mmap_mutex = PTHREAD_MUTEX_INITIALIZER;

/* returns 0 if an empty interval was found else the interval size,
   lindex and rindex are only set if return value is positive */
int IntervalQuery::find_interval(unsigned *pos_map, off_t num_entries,
		   unsigned begin, unsigned end, unsigned *lindex, unsigned *rindex, unsigned left_limit, unsigned right_limit)
{
  //printf("find_interval: %d %d\n",begin,end);

	unsigned val, mid, left = left_limit, right = right_limit;

	//assert(begin <= end);
	if (begin>end)
	{
		fprintf(stderr, "[find_interval] Warning: begin=%i > end=%i\n", begin, end) ;
		return 0 ;
	}

	if (num_entries==0)
	{
		*lindex=0 ;
		*rindex=0 ;
		return 0 ;
	}
	
	if (pos_map[right] < begin || pos_map[left] > end) {
#ifdef DEBUG
	  printf("ret 0\n\n");
#endif
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
	//left = left_limit;
	right = right_limit;
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
        //assert(*rindex >= *lindex);
        return *rindex - *lindex + 1;
}


int IntervalQuery::mmap_file(const char *filename, int open_mode, void **map, off_t *size, bool convert)
{
	static size_t mmap_total_size = 0 ;
	
	std::string filename_s ;

	filename_s.assign(filename) ;
	pthread_mutex_lock(&mmap_mutex) ;
	
	for (size_t i=0; i<mmap_fname_list.size(); i++)
		if (filename_s==mmap_fname_list[i])
		{
			//fprintf(stdout, "found item %i (%s)\n", (int) i, filename) ;
			
			*map=mmap_ptr_list[i] ;
			*size=mmap_size_list[i] ;

			pthread_mutex_unlock(&mmap_mutex) ;

			return 1 ;
		}
	pthread_mutex_unlock(&mmap_mutex) ;
	

	int fd, ret, mmap_prot, mmap_flags;
	struct stat file_status;

	if (open_mode == O_RDONLY) {
		mmap_prot = PROT_READ;
		mmap_flags = MAP_PRIVATE;
	} else {
		mmap_prot = PROT_READ | PROT_WRITE;
		mmap_flags = MAP_SHARED;
	}
	/*printf("filename: %s\n", filename);*/
	ret = open(filename, open_mode);
	if (ret < 0)
	{
		char buf[1000] ;
		sprintf(buf, "can not open %s: %s\n", filename, strerror(errno));
		throw IntervalQueryException(buf) ;
	}

	fd = ret;
	if (fstat(fd, &file_status) < 0)
	{
		char buf[1000] ;
		sprintf(buf, "fstat error: %s\n", strerror(errno)) ;
		throw IntervalQueryException(buf) ;
	}
	*size = file_status.st_size;
	if (!*size)
	{
		fprintf(stdout, "mmap warning: %s is empty\n", filename) ;
		*map = NULL ;
		ret = 1 ;
		close(fd) ;
		return ret ;
		//throw IntervalQueryException(buf) ;
	}
	*map = mmap(NULL, *size, mmap_prot, mmap_flags, fd, 0);
	//fprintf(stderr, "##mmap %ld size %ld\n", (size_t)*map, *size) ;
	mmap_total_size+=*size ;
	
	if (*map == MAP_FAILED)
	{
		fprintf(stderr, "mmap error: %s\n", strerror(errno));
		close(fd);
	    return -1 ;
	}
	ret = 1;
	close(fd);

	//fprintf(stdout, "adding mmap #%i (%s)\n", (int)mmap_fname_list.size(), filename) ;

	//fprintf(stderr, "mmap total size: %ld \n", mmap_total_size) ;
	if (convert)
	{
		// convert float score files to short, unmap memory and use new memory block
		float *map_f = (float*) *map ;
		unsigned short int *map_i = (unsigned short int*) malloc((*size*sizeof(unsigned short int))/sizeof(float)) ;
		for (size_t i=0; i< *size/sizeof(float); i++)
		{
			assert(map_f[i]>=0.0 && map_f[i]<=1.0) ;
			map_i[i] = map_f[i]*65535 ;
		}
		munmap(*map, *size);
		//fprintf(stderr, "##munmap %ld size %ld\n", (size_t)*map, *size) ;
		mmap_total_size = mmap_total_size - *size ;
		//fprintf(stderr, "mmap total size: %ld -\n", mmap_total_size) ;

		*size = ((*size)*sizeof(unsigned short int))/sizeof(float) ;
		*map = map_i ;

		mmap_total_size+=*size ;
		//fprintf(stderr, "mmap total size: %ld +\n", mmap_total_size) ;
	}

	pthread_mutex_lock(&mmap_mutex) ;
	mmap_fname_list.push_back(filename_s) ;
	mmap_ptr_list.push_back(*map) ;
	mmap_size_list.push_back(*size) ;
	pthread_mutex_unlock(&mmap_mutex) ;
	
	return ret;
}


/*
 * rhs[0] basename
 * rhs[1] list of m score names
 * rsh[2] Interval Matrix (n x 2)
 *
 * return:
 * lhs[0]: list of positions that lie in any given interval, k := number of such positions
 * lhs[1]: k x m score Matrix
 * lhs[2]: index vector (k x 1), lhs[2][j] = l <=> rhs[2][l, 0] <= lhs[0][j] <= rhs[2][l, 1].
 */

int IntervalQuery::interval_query(char* basename, char** score_names, unsigned int num_scores, int* interval_matrix, unsigned int num_intervals)
{
	clock_t start_time=clock() ;

   //printf("Entering interval_query...\n");
	char filename[MAXLINE], *err_txt;
	int ret;
	unsigned i, j = 0;

#ifdef CONVERT_SCORE_MAPS
	unsigned short int **score_maps = NULL;
	assert(sizeof(unsigned short int)==2) ;
#else
	float **score_maps = NULL;
#endif
	
	off_t *score_sizes = NULL;
	//double *interval_ptr;
	int* interval_ptr;
	unsigned m, k,total_num_pos, *pos_map = NULL;
	off_t pos_size, num_positions;
	unsigned llimit=0;
	unsigned rlimit=0;
	unsigned left;
	unsigned right;
	unsigned lindex, rindex;


   //printf("initialized local vars...\n");

	ret = -E_SYNTAX;
	err_txt = (char*)"fatal: no interval(s) given";
	if (!num_intervals)
		goto out;
	err_txt = (char*)"fatal: no score file(s) given";
	if (!num_scores)
		goto out;

#ifdef CONVERT_SCORE_MAPS
	score_maps = new unsigned short int*[num_scores] ;//(float**) malloc(num_scores * sizeof(void *));
#else
	score_maps = new float*[num_scores] ;//(float**) malloc(num_scores * sizeof(void *));
#endif
	score_sizes = new off_t[num_scores] ;//(off_t*) malloc(num_scores * sizeof(off_t *));

	memset(score_maps, 0, num_scores * sizeof(void *));

   //printf("setting score_maps...\n");

	//printf("%s.pos\n", basename);
	sprintf(filename, "%s.pos", basename);
	ret = mmap_file(filename, O_RDONLY, (void **)&pos_map, &pos_size, false);
	err_txt = (char*)"pos mmap error";
	if (ret < 0)
		goto out;
	num_positions =	pos_size / sizeof(unsigned);

   //printf("Loading file / mmapping...\n");

#ifdef DEBUG
	err_txt = (char*)"pos map is not sorted";
	printf("num_positions: %i\n", num_positions);
	for (i=1; i<num_positions; ++i) {
	    if (pos_map[i-1] > pos_map[i])
	            goto out;
	}
	printf("successfully checked that pos_map is sorted\n");
#endif

	for (i = 0; i < num_scores; i++) {
		//printf("%s.%s\n", basename, score_names[i]);
		sprintf(filename, "%s.%s", basename, score_names[i]);
#ifdef CONVERT_SCORE_MAPS
		ret = mmap_file(filename, O_RDONLY, (void **)&score_maps[i], &score_sizes[i], true);
#else
		ret = mmap_file(filename, O_RDONLY, (void **)&score_maps[i], &score_sizes[i], false);
#endif
		err_txt = (char*)"score mmap error";
		if (ret < 0)
		        goto out;
	}
	total_num_pos = 0;

   interval_ptr = interval_matrix;
   assert(interval_ptr[0]>=0) ;
   left = (unsigned)interval_ptr[0];
   assert(interval_ptr[2 * (num_intervals-1) + 1]>=0) ;
   right = (unsigned)interval_ptr[2 * (num_intervals-1) + 1];
   ret = find_interval(pos_map, num_positions, left, right, &lindex, &rindex, 0, num_positions-1);
   if (ret>0)
   {
	   llimit=lindex;
	   rlimit=rindex;
   }

	for (j = 0; j < num_intervals; j++) {
	    left = (unsigned)interval_ptr[2 * j];
	    right = (unsigned)interval_ptr[2 * j + 1];
       //printf("left/right are %d/%d\n",left,right);

#ifdef DEBUG
	      for (i=0; i<num_positions; ++i) {
		if (left <= pos_map[i] && pos_map[i] <= right)
		  printf("%i ", pos_map[i]);
	      }
	      printf("\n\n");
#endif

	    ret = -E_SYNTAX;
	    err_txt = (char*)"bad interval";
	    if (left > right)
	        goto out;
		if (num_positions>1)
			ret = find_interval(pos_map, num_positions, left, right, &lindex, &rindex, llimit, rlimit);
	    if (ret > 0) {/* non-empty interval */
	        total_num_pos += rindex - lindex + 1;
#ifdef DEBUG
		printf("lidx: %u, ridx: %u, t: %u\n", lindex, rindex, total_num_pos);
#endif
	    }
	}

   //printf("after 1st call to find_interval\n");

   //pos_ptr = (int*) malloc(total_num_pos*sizeof(int));
	//score_ptr = (double*) malloc(total_num_pos*num_scores*sizeof(double));
   //index_ptr = (int*) malloc(total_num_pos*sizeof(int));

   pos_ptr = new int[total_num_pos];
   score_ptr = new double[total_num_pos*num_scores];
   index_ptr = new int[total_num_pos];

   num_hits_p = total_num_pos;

	k = 0;
	for (j = 0; j < num_intervals; j++) 
	{
	    unsigned left_ = (unsigned)interval_ptr[2 * j];
	    unsigned right_ = (unsigned)interval_ptr[2 * j + 1];
	    unsigned lindex_, rindex_;
		
	    ret = find_interval(pos_map, num_positions, left_, right_, &lindex_, &rindex_, llimit, rlimit);
	    if (ret <= 0) /* empty interval */
		continue;
		assert((unsigned)ret == rindex_ - lindex_ + 1);
	    for (m = 0; m < (unsigned)ret; m++) {
	      if (index_ptr)
	  	index_ptr[k + m] = j;
	      err_txt = (char*)"returned position out of interval";
              if (left_ > pos_map[lindex_ + m] || pos_map[lindex_ + m] > right_)
		goto out;

#ifdef DEBUG
              printf("l: %u, r: %u, p: %u\n", left_, right_, pos_map[lindex_ + m]);
#endif
	      pos_ptr[k + m] = pos_map[lindex_ + m];
	      for (i = 0; i < num_scores; i++)
		  {
			  /*score_ptr[i + num_scores * (k + m)] = (double)score_maps[i][lindex_ + m];*/
		    if ((lindex_ + m)*sizeof(score_maps[i][0]) >= (unsigned) score_sizes[i])
		      score_ptr[i * total_num_pos + (k + m)] = 0.0 ; // rescue, when the conf_cum file is too short
		    else
		      {
#ifdef CONVERT_SCORE_MAPS
			score_ptr[i * total_num_pos + (k + m)] = ((double) score_maps[i][lindex_ + m])/65535 ;
#else
			score_ptr[i * total_num_pos + (k + m)] = (double) score_maps[i][lindex_ + m] ;
#endif
		      }
		  }
		}
	    k += ret;
	}

	ret = 1;

   //printf("after 2nd call to find_interval\n");

 out:

	pthread_mutex_lock(&mmap_mutex) ;
	total_time += clock()-start_time ;
	pthread_mutex_unlock(&mmap_mutex) ;
	
	/*for (i = 0; i < num_scores; i++) {
	  if (score_maps[i])
	    munmap(score_maps[i], score_sizes[i]);
		}*/

	delete[] score_maps ;
	delete[] score_sizes ;
	/*if (pos_map)
	  munmap(pos_map, pos_size);*/

	return ret ;

   //printf("Leaving interval_query...\n");
}


void IntervalQuery::cleanup()
{
         if (pos_ptr != NULL)
		 {
			 delete[] pos_ptr;
			 pos_ptr=NULL ;
		 }

         if (index_ptr != NULL)
		 {
			 delete[] index_ptr;
			 index_ptr = NULL ;
		 }

         if (score_ptr != NULL)
		 {
			 delete[] score_ptr;
			 score_ptr=NULL ;
		 }

}

int IntervalQuery::query(char* basename, char** score_names, unsigned int num_scores, int* interval_matrix, unsigned int num_intervals)
{
	int ret = interval_query(basename, score_names, num_scores, interval_matrix, num_intervals);
	if (ret<0)
		return ret ;

	num_scores_p = num_scores;
	return num_hits_p;
}


void IntervalQuery::getResults(int* pos, int* index, double* score) 
{

   for(unsigned int k=0; k<num_hits_p;k++) {
      pos[k] = pos_ptr[k];
      index[k] = index_ptr[k];
   }

   for(unsigned int k=0; k<num_hits_p*num_scores_p;k++) {
      score[k] = score_ptr[k];
   }

}
