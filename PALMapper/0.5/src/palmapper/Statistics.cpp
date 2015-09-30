#include <palmapper/Statistics.h>
#include <assert.h>

Statistics::Statistics() 
{
	PERFECT_READS = 0;
	PERFECT_HITS = 0;
	PERFECT_HITS_REV = 0;
	NUM_HITS = 0;
	for (unsigned int i = 0; i < Config::MAX_READ_LENGTH; ++i)
		HITS_LEN[i] = 0;
	for (int i = 0; i != Config::MAX_EDIT_OPS + 1; ++i)
		HITS_MM[i] = 0;
	READS_MAPPED = 0;
	NUM_ALIGNMENTS = 0;
	NUM_WHOLE_ALIGNMENTS = 0;
	ENDSTART_MAPPED[0] = 0;
	ENDSTART_MAPPED[1] = 0;
	NOT_ALIGNED[0] = 0;
	NOT_ALIGNED[1] = 0;
	NUM_READS = 0;
	HITS_PER_READ = 0;
	GAPS_ENCOUNTERED[0] = 0;
	GAPS_ENCOUNTERED[1] = 0;
	GAPS_ENCOUNTERED[2] = 0;
	TOO_MANY_MMS[0] = 0;
	TOO_MANY_MMS[1] = 0;
	BREAK_GLOBAL_ALIGNMENT[0] = 0;
	BREAK_GLOBAL_ALIGNMENT[1] = 0;
	BREAK_TB_IN_GLOBAL_ALIGNMENT = 0;
	CELLS_GLOBAL = 0;
	CELLS_OVERHANG = 0;
	W = 0;
	listcount = 0;
	listocc = 0;

	//Alignments
	int ret = pthread_mutex_init(&alignment_num_mutex, NULL) ;
	assert(ret==0) ;
	alignment_num_unmapped=0 ;
	alignment_num_unspliced_best=0 ;
	alignment_num_unspliced_suboptimal=0 ;
	alignment_num_spliced_best=0 ;
	alignment_num_spliced_suboptimal=0 ;
	alignment_last_spliced_report = clock() ;

	//QPALMA Filter
	qpalma_last_filter_report=0;
	for (int i=0; i<num_filter_reasons; i++)
	{
		qpalma_filter_stat_spliced[i]=0 ;
		qpalma_filter_stat_unspliced[i]=0 ;
	}
	

	//QPALMA
	qpalma_last_timing_report=0 ;
	qpalma_region1_time=0; 
	qpalma_region_align_time=0 ;
	qpalma_align_time=0 ;
	intervalquery_total_time=0 ;
	
	qpalma_total_dna_length=0 ;
	qpalma_total_alignments=0 ;
	qpalma_total_num_threads=0 ;
	qpalma_total_num_thread_tasks=0 ;


	// HITS
	hits_seed2genome=0, hits_last_report_total=0;
	hits_seed2genome_cnt=0;

	hits_seek=0, hits_part1=0, hits_part2=0, hits_part3=0, hits_part4=0, hits_part5=0,hits_last_report=0 ;
	hits_seek_cnt=0, hits_part1_cnt=0, hits_part2_cnt=0, hits_part3_cnt=0, hits_part4_cnt=0, hits_part5_cnt=0 ;

	ret = pthread_mutex_init(&hit_seed_mutex, NULL) ;
	assert(ret==0) ;

}
