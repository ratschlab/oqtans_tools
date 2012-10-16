/*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 3 of the License, or
*  (at your option) any later version.
*
*   Written (W) 2010-2011 Jonas Behr, Regina Bohnert, Gunnar Raetsch
*   Copyright (C) 2010-2011 Max Planck Society
*/


#ifndef __READ_H__
#define __READ_H__

#include <stdint.h>
#include <cctype>
#include <stdio.h>
#include <assert.h>
#include <vector>
		using std::vector;


class CRead {		 
 public:
		/** constructor
		 */
		CRead();
		~CRead();
		
		vector<int> block_starts;
		vector<int> block_lengths;
		char* read_id;
		char* sam_line;
		int start_pos;
		char * strand;
		int matches;
		int mismatches;
		int multiple_alignment_index;
		bool left; 
		bool right;
		bool reverse;
		bool is_clipped;
		
		void get_coverage(int p_start_pos, int p_end_pos, uint32_t* coverage);
		int get_last_position();
		void get_reads_sparse(int p_start_pos, int p_end_pos, double* reads, uint32_t & reads_c, uint32_t row_idx);
		void get_introns(vector<int>* introns);
		void get_introns(vector<uint32_t>* intron_starts, vector<uint32_t>* intron_ends, vector<uint32_t>* block_len1, vector<uint32_t>* block_len2); 
		void get_acc_splice_sites(vector<int>* acc_pos);
		void get_don_splice_sites(vector<int>* acc_pos);
		int max_intron_len();
		int min_exon_len();
		bool operator==(const CRead& read) const;
		void print();
		void set_strand(char s);
		int get_mismatches();

		bool is_spliced()
		{
			return block_starts.size()>1;
		};

		/** check if this read is compatible with a given exon */
		bool is_compatible(int exon_start, int exon_end)
		{	
			assert(exon_start<exon_end);
		
			// case: no overlapp 1
			if (exon_end<start_pos)
				return false;
			// case: no overlapp 2
			if (exon_start>get_last_position())
				return false;

			// case: mapped read has to be completely within the exon 
			if (block_starts.size()==1)
				return (exon_start<=start_pos) && (exon_end>=get_last_position()); 

			for (uint32_t i=0; i<block_starts.size(); i++)
			{
				int start = block_starts[i]+start_pos;
				int end = block_starts[i]+start_pos+block_lengths[i];

				bool first = i==0;
				bool last = i==block_starts.size()-1;

				if (first && start>=exon_start && end==exon_end)
					return true;
				if (last && start == exon_start && end<=exon_end)
					return true;
				// case two introns in one read
				if (start == exon_start && end == exon_end)
					return true;
			}
			return false;
		};

		static bool compare_by_read_id(const CRead* read1, const CRead* read2)
		{
			if (!read1->read_id)
				return true;
			if (!read2->read_id)
				return false;
		
			int cnt1=0;
			while (read1->read_id[cnt1]!='\0')
				cnt1++;
			int cnt2 = 0;
			while (read2->read_id[cnt2]!='\0')
				cnt2++;
			
			return std::lexicographical_compare(read1->read_id,read1->read_id+cnt1,read2->read_id,read2->read_id+cnt2);
		};
};
#endif
