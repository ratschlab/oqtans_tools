/*

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Written (W) 2007-2009 Fabio De Bona, Andre Noll, Gunnar Raetsch
Copyright (C) 2007-2009 Max-Planck-Society

*/

#ifndef __INTERVALQUERY_H__
#define __INTERVALQUERY_H__

#include <vector>
#include <string>

#include <stdexcept> // stdexcept header file contains runtime_error
using std::runtime_error;

class IntervalQueryException : public runtime_error
{
public:
	IntervalQueryException(const char * errmsg)
		: runtime_error(errmsg) 
		{
		}
} ;

class IntervalQuery {

   private:
      int *pos_ptr;
      int *index_ptr;
      double *score_ptr;

      unsigned int num_scores_p;
      unsigned int num_hits_p;

	  static std::vector<std::string> mmap_fname_list ;
	  static std::vector<void*> mmap_ptr_list ;
	  static std::vector<int> mmap_size_list ;
	  static pthread_mutex_t mmap_mutex ;

public:	  
	  static clock_t total_time ;
	  
   public:
      IntervalQuery() {
         pos_ptr = 0;
         index_ptr = 0;
         score_ptr = 0;
         num_scores_p  = 0;
         num_hits_p    = 0;
      }
      
      ~IntervalQuery() {
		  cleanup() ;
	  }
   
	  void cleanup() ;
	  
	  int find_interval(unsigned *pos_map, off_t num_entries, unsigned begin, unsigned end, unsigned *lindex, unsigned *rindex, unsigned left_limit, unsigned right_limit) ;
	  
	  int mmap_file(const char *filename, int open_mode, void **map, off_t *size, bool convert) ;
	  
      int interval_query(char* basename, char** score_names, unsigned int num_scores, int* interval_matrix, unsigned int num_intervals);

      int query(char* basename, char** score_names, unsigned int num_scores, int* interval_matrix, unsigned int num_intervals);

      void getResults(int* pos, int* index, double* score);

	  int getResultSize()
		  {
			  return num_hits_p ;
		  }
};

#endif // __INTERVALQUERY_H_
