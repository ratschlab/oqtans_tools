// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "palmapper.h"
#include "print.h"
#include <assert.h>
#include <vector>

#include <palmapper/Genome.h>
#include <palmapper/Hits.h>
#include <palmapper/QueryFile.h>
#include <palmapper/print.h>

static const int MAX_EXON_LEN=100 ;

void print_stats(QueryFile &queryFile) {
	printf("\n########### _config.STATISTICS ###########\n");
	printf("Mapped Reads: %i of all %d reads\n", _stats.READS_MAPPED, _stats.NUM_READS);
	unsigned int i;
	for (i = 0; i != (unsigned int)_config.NUM_MISMATCHES + 1; ++i)
		printf(" Reads with %d mismatches: %u\n", i, _stats.HITS_MM[i]);
	printf(
			"  Perfect Plus-Hits: %i\n  Perfect Minus-Hits:\t%i\n  Total:\t\t%i\n",
			_stats.PERFECT_HITS, _stats.PERFECT_HITS_REV, _stats.PERFECT_HITS + _stats.PERFECT_HITS_REV);
	printf("   Perfect matching reads (+ or -): %i\n", _stats.PERFECT_READS);

	for (i = _config.INDEX_DEPTH; i != queryFile.max_length(); ++i)
		printf("    Total Hits of length %d: %u - %.4f%%\n", i, _stats.HITS_LEN[i],
				((double) (100 * _stats.HITS_LEN[i]) / _stats.NUM_HITS));
	printf("\n    Total number of hits:                  %lu\n", _stats.NUM_HITS);
	printf("      Number of overhang alignments:       %lu\n", _stats.NUM_ALIGNMENTS);
	printf("      Successful overhang alignments:      %lu\n", _stats.NUM_ALIGNMENTS
			- _stats.GAPS_ENCOUNTERED[0] - _stats.GAPS_ENCOUNTERED[1] - _stats.GAPS_ENCOUNTERED[2]
			- _stats.TOO_MANY_MMS[0] - _stats.TOO_MANY_MMS[1]);
	printf("      Number of global alignments:         %lu\n",
			_stats.NUM_WHOLE_ALIGNMENTS);
	printf("      Successful global alignments:        %lu\n",
			_stats.NUM_WHOLE_ALIGNMENTS - _stats.BREAK_GLOBAL_ALIGNMENT[0]
					- _stats.BREAK_TB_IN_GLOBAL_ALIGNMENT);
	printf("      Breaks in global alignments(no hit): %lu\n",
			_stats.BREAK_GLOBAL_ALIGNMENT[0]);
	printf("      Breaks in global alignments(hit):    %lu\n",
			_stats.BREAK_GLOBAL_ALIGNMENT[1]);
	printf("      Breaks after global alignment:       %lu\n\n",
			_stats.BREAK_TB_IN_GLOBAL_ALIGNMENT);
	printf("      Hits of len %d, not aligned:         %d\n", queryFile.max_length() - 1, _stats.NOT_ALIGNED[0]);
	printf("      Mappings of fast mapping:            %d\n", _stats.NOT_ALIGNED[1]);
	//printf("      hits of len 36, not aligned (with 1MM at start/beg of read): %d\n",_stats.NOT_ALIGNED[1]);
	//printf("        Number of unmapped hits at beg/end of chrom: %d and in align-step: %d\n", _stats.ENDSTART_MAPPED[0], _stats.ENDSTART_MAPPED[1]);
	printf("      Overlapping hits at beg/end of chr:  %d\n",
			_stats.ENDSTART_MAPPED[0]);

	printf("\nGaps encountered on diagonal during overhang-alignment:   %lu\n",
			_stats.GAPS_ENCOUNTERED[0]);
	printf("Gaps encountered due to score during overhang-alignment:  %lu\n",
			_stats.GAPS_ENCOUNTERED[1]);
	printf("Gaps encountered after overhang-alignment:                %lu\n",
			_stats.GAPS_ENCOUNTERED[2]);
	printf("Too many MMs during overhang-alignment (discarded hits):  %lu\n",
			_stats.TOO_MANY_MMS[0]);
	printf("Too many MMs after overhang-alignment (discarded hits):   %lu\n\n",
			_stats.TOO_MANY_MMS[1]);

	printf("cells processed in overhang: %lu\n", _stats.CELLS_OVERHANG);
	printf("cells processed in global:   %lu\n", _stats.CELLS_GLOBAL);
	printf("total cells:                 %lu\n\n", _stats.CELLS_OVERHANG
			+ _stats.CELLS_GLOBAL);

	printf("gap limit overcome: %lu\n\n", _stats.W);

	//printf("reads filtered out due to too much non-base chars: %i\n\n",READS_FILTERED);

	/*printf("Print in stat... ");
	 unsigned int sum=0;
	 FILE* f;
	 f = fopen("stat","w");
	 for (i=0; i!= _config.MAX_READ_LENGTH-_config.INDEX_DEPTH+1; ++i) {
	 for (j=0; j!=_config.MAX_READ_LENGTH - _config.INDEX_DEPTH - i + 1; ++j) {
	 if (*(*(_stats.HITS_READPOS+i)+j) != 0) {
	 fprintf(f, "\t(%d)%u ", j+1, *(*(_stats.HITS_READPOS+i)+j));
	 sum += *(*(_stats.HITS_READPOS+i)+j);
	 }
	 }
	 fprintf(f, "\nLENGTH %d: %d", i+_config.INDEX_DEPTH, sum);
	 if (i!=_config.MAX_READ_LENGTH-_config.INDEX_DEPTH && sum != 0 && i<36) fprintf(f, "\tfirst: %.4f%%, last: %.4f%%, other: %.4f%%\n", ((double) (100 * *(*(_stats.HITS_READPOS+i))) / sum), ((double) (100 * (*(*(_stats.HITS_READPOS+i)+(36 - _config.INDEX_DEPTH - i))))/sum), 100 - ((double) (100 * *(*(_stats.HITS_READPOS+i))) / sum) - ((double) (100 * (*(*(_stats.HITS_READPOS+i)+(36 - _config.INDEX_DEPTH - i))))/sum));
	 else fprintf(f, "\n");
	 sum = 0;
	 }
	 printf("done\n");*/
}

/*int different_hit(HIT* hit1, HIT* hit2)
{
	//fprintf(stdout, "diff %i %i\n", hit1->start, hit2->start) ;
	
	if (abs(hit1->start-hit2->start)<(int)READ_LENGTH)
		return 0 ;
	return 1 ;
}*/




// sorts by position, if it's identical, then gaps have to be before mismatches
int compare_editops(const void *a, const void *b) {
	EDIT_OPS *ia = (EDIT_OPS *) a;
	EDIT_OPS *ib = (EDIT_OPS *) b;
	if (ia->pos == ib->pos)
		return (int) (ib->mm - ia->mm);

	return (int) (abs(ia->pos) - abs(ib->pos));
}

int compare_int(const void *a, const void *b) 
{
	int *ia = (int *) a;
	int *ib = (int *) b;

	return (int) (*ia - *ib);
}

void print_leftovers(Read const &read, const char * tag, std::ostream *LEFTOVER_FP)
{
	_stats.alignment_num_unmapped++ ;
	
	if (read.get_orig()==NULL)
	{
		if (read.format() == 0)
			fprintf(LEFTOVER_FP, "@%s%s\n%s\n+\n%s\n", read.id(), tag, read.data(), read.quality(0));
		else if (read.format() == 1)
			fprintf(LEFTOVER_FP, ">%s%s\n%s\n", read.id(), tag, read.data());
		else {
			fprintf(LEFTOVER_FP, "%s%s\t%s\t%d", read.id(), tag, read.data(), read.pe_flag());
			if (strlen(read.quality(0)) > 0) fprintf(LEFTOVER_FP, "\t%s", read.quality(0));
			if (strlen(read.quality(1)) > 0) fprintf(LEFTOVER_FP, "\t%s", read.quality(1));
			if (strlen(read.quality(2)) > 0) fprintf(LEFTOVER_FP, "\t%s", read.quality(2));
			fprintf(LEFTOVER_FP, "\n");
		}
	}
	else
	{
		if (read.format() == 0)
			fprintf(LEFTOVER_FP, "@%s%s\n%s\n+\n%s\n", read.id(), tag, read.get_orig()->data(), read.get_orig()->quality(0));
		else if (read.format() == 1)
			fprintf(LEFTOVER_FP, ">%s%s\n%s\n", read.id(), tag, read.get_orig()->data());
		else {
			fprintf(LEFTOVER_FP, "%s%s\t%s\t%d", read.id(), tag, read.get_orig()->data(), read.pe_flag());
			if (strlen(read.get_orig()->quality(0)) > 0) fprintf(LEFTOVER_FP, "\t%s", read.get_orig()->quality(0));
			if (strlen(read.get_orig()->quality(1)) > 0) fprintf(LEFTOVER_FP, "\t%s", read.get_orig()->quality(1));
			if (strlen(read.get_orig()->quality(2)) > 0) fprintf(LEFTOVER_FP, "\t%s", read.get_orig()->quality(2));
			fprintf(LEFTOVER_FP, "\n");
		}
	}
}

void print_alignment_matrix(Read const &read, int chrstart, int readstart, int length,
		int offset_front, int offset_end, Chromosome const &chr, char ori, int K) {
	int i, j;

	printf(" k=%d |\t-\t", K);
	if (ori == '+')
		for (i = 0; i != length; i++)
			printf("%c\t", chr[chrstart + i]);
	else
		for (i = 0; i != length; i++)
			printf("%c\t", get_compl_base(chr[chrstart - i]));
	printf(
			"\n----------------------------------------------------------------------------------------------");
	printf(
			"----------------------------------------------------------------------------------------------\n");

	for (j = 0; j != length + 1; ++j) {
		if (j == 0)
			printf("  -  |\t");
		else {
			if ((readstart == 0 && (j <= offset_front || (j > length
					- offset_end && K != 1))) || (readstart != 0 && j > length
					- offset_end))
				printf("  X  |\t");
			else
				printf("  %c  |\t", read.data()[readstart + j - (readstart == 0)
						* offset_front - 1]);
		}

		if (j > K)
			for (i = 0; i != j - K; ++i)
				printf("\t");

		for (i = 0; i != 2* K + 1; ++i) {
			if (i - j > K)
				break;
			if (i + j > length + K)
				break;
			if (i > length || j > length)
				break;
			if ((readstart != 0 && j > length - offset_end) || (readstart == 0
					&& j < offset_front))
				break;
			//printf("{i%d,j%d}",i,j);
			//printf("%.1f(%c)\t", M[i][j], T[i][j]);
			if (i == 6 && j == 1)
				printf("---");
		}

		printf(
				"\n----------------------------------------------------------------------------------------------");
		printf(
				"----------------------------------------------------------------------------------------------\n");
	}

}

