#include <stdio.h> 
#include <unistd.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "bwtaln.h"
#include "bwtgap.h"
#include "utils.h"
#include "bntseq.h"
#include "bwtmyaln.h"
#include <assert.h>

gap_opt_t *mygap_init_opt()
{
	gap_opt_t *o;
	o = (gap_opt_t*)calloc(1, sizeof(gap_opt_t));
	/* IMPORTANT: s_mm*10 should be about the average base error
	   rate. Voilating this requirement will break pairing! */
	o->s_mm = 3; o->s_gapo = 11; o->s_gape = 4;
	o->max_diff = 0; o->max_gapo = 0; o->max_gape = 0;
	o->indel_end_skip = 5; o->max_del_occ = 10; o->max_entries = 2000000;
	o->mode = BWA_MODE_GAPE | BWA_MODE_COMPREAD;
	o->seed_len = 32; o->max_seed_diff = 2;
	o->fnr = 0.0;
	o->n_threads = 1;
	o->max_top2 = 30;
	o->trim_qual = 0;
	return o;
}

int bwa_cal_maxdiff(int l, double err, double thres) ;

static int bwt_cal_width(const bwt_t *rbwt, int len, const ubyte_t *str, bwt_width_t *width)
{
	bwtint_t k, l, ok, ol;
	int i, bid;
	bid = 0;
	k = 0; l = rbwt->seq_len;
	for (i = 0; i < len; ++i) {
		ubyte_t c = str[i];
		if (c < 4) {
			bwt_2occ(rbwt, k - 1, l, c, &ok, &ol);
			//fprintf(stdout, "i=%i, c=%i, k=%i, l=%i, ok=%i, ol=%i\n", i, c, k, l, ok, ol) ;
			k = rbwt->L2[c] + ok + 1;
			l = rbwt->L2[c] + ol;
		}
		if (k > l || c > 3) { // then restart
			k = 0;
			l = rbwt->seq_len;
			++bid;
		}
		width[i].w = l - k + 1;
		width[i].bid = bid;
	}
	width[len].w = 0;
	width[len].bid = ++bid;
	return bid;
}

void mybwa_cal_sa_reg_gap(int tid, bwt_t *const bwt[2], int n_seqs, bwa_seq_t *seqs, const gap_opt_t *opt)
{
	int i, max_l = 0, max_len;
	gap_stack_t *stack;
	bwt_width_t *w[2], *seed_w[2];
	const ubyte_t *seq[2];
	gap_opt_t local_opt = *opt;

	// no mismatches or gaps in seed
	local_opt.max_gape=0 ;
	local_opt.max_gapo=0 ;
	local_opt.max_seed_diff=0 ;
	local_opt.max_diff=0 ;

	// initiate priority stack
	for (i = max_len = 0; i != n_seqs; ++i)
		if (seqs[i].len > max_len) max_len = seqs[i].len;
	if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(max_len, BWA_AVG_ERR, opt->fnr);
	if (local_opt.max_diff < local_opt.max_gapo) local_opt.max_gapo = local_opt.max_diff;
	stack = gap_init_stack(local_opt.max_diff, local_opt.max_gapo, local_opt.max_gape, &local_opt);

	seed_w[0] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	seed_w[1] = (bwt_width_t*)calloc(opt->seed_len+1, sizeof(bwt_width_t));
	w[0] = w[1] = 0;
	for (i = 0; i != n_seqs; ++i) {
		bwa_seq_t *p = seqs + i;
		p->sa = 0; p->type = BWA_TYPE_NO_MATCH; p->c1 = p->c2 = 0; p->n_aln = 0; p->aln = 0;
		seq[0] = p->seq; seq[1] = p->rseq;
		if (max_l < p->len) {
			max_l = p->len;
			w[0] = (bwt_width_t*)realloc(w[0], (max_l + 1) * sizeof(bwt_width_t));
			w[1] = (bwt_width_t*)realloc(w[1], (max_l + 1) * sizeof(bwt_width_t));
			memset(w[0], 0, (max_l + 1) * sizeof(bwt_width_t));
			memset(w[1], 0, (max_l + 1) * sizeof(bwt_width_t));
		}
		bwt_cal_width(bwt[0], p->len, seq[0], w[0]);
		bwt_cal_width(bwt[1], p->len, seq[1], w[1]);
		//fprintf(stdout, "w[0]=%i, w[1]=%i\n", w[0]->w, w[1]->w) ;
		
		if (opt->fnr > 0.0) local_opt.max_diff = bwa_cal_maxdiff(p->len, BWA_AVG_ERR, opt->fnr);
		local_opt.seed_len = opt->seed_len < p->len? opt->seed_len : 0x7fffffff;
		if (p->len > opt->seed_len) {
			bwt_cal_width(bwt[0], opt->seed_len, seq[0] + (p->len - opt->seed_len), seed_w[0]);
			bwt_cal_width(bwt[1], opt->seed_len, seq[1] + (p->len - opt->seed_len), seed_w[1]);
		}
		// core function
		p->aln = bwt_match_gap(bwt, p->len, seq, w, p->len <= opt->seed_len? 0 : seed_w, &local_opt, &p->n_aln, stack);
		// store the alignment
		
		free(p->name); free(p->seq); free(p->rseq); free(p->qual);
		p->name = 0; p->seq = p->rseq = p->qual = 0;
	}
	
	free(seed_w[0]); free(seed_w[1]);
	free(w[0]); free(w[1]);
	gap_destroy_stack(stack);
}


void bwa_sai2sam_se_core(const char *prefix, const char *fn_sa, const char *fn_fa, int n_occ) ;

int bns_coor_pac2real(const bntseq_t *bns, int64_t pac_coor, int len, int32_t *real_seq) ;
extern "C" int64_t pos_end(const bwa_seq_t *p) ;


void mybwa_cal_pac_pos_core(const bwt_t *forward_bwt, const bwt_t *reverse_bwt, bwa_seq_t *seq, const int max_mm, const float fnr)
{
	int max_diff;
	if (seq->type != BWA_TYPE_UNIQUE && seq->type != BWA_TYPE_REPEAT) return;
	max_diff = fnr > 0.0? bwa_cal_maxdiff(seq->len, BWA_AVG_ERR, fnr) : max_mm;
	if (seq->strand) { // reverse strand only
		//fprintf(stdout, "forward used\n") ;
		seq->pos = bwt_sa(forward_bwt, seq->sa);
		//seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
	} else { // forward strand only
		/* NB: For gapped alignment, p->pos may not be correct, which
		 *     will be fixed in refine_gapped_core(). This line also
		 *     determines the way "x" is calculated in
		 *     refine_gapped_core() when (ext < 0 && is_end == 0). */
		//fprintf(stdout, "backward used\n") ;
		seq->pos = reverse_bwt->seq_len - (bwt_sa(reverse_bwt, seq->sa) + seq->len);
		//seq->pos = bwt_sa(reverse_bwt, seq->sa);
		//seq->seQ = seq->mapQ = bwa_approx_mapQ(seq, max_diff);
	}

	//fprintf(stdout, "seq->pos=%i\n", seq->pos) ;
}


static const gap_opt_t* bwt_opt=NULL ;
static bwt_t *bwt_bwt[2]={NULL, NULL};
static bntseq_t *bwt_bns = NULL;

extern "C" void bwa_seed2genome_init(const char *prefix, gap_opt_t *opt)
{
	if (!opt)
		opt=mygap_init_opt() ;
	opt->mode=BWA_MODE_BAM_SE ;
	bwt_opt = opt ;
	
	{ // load BWT
		char *str = (char*)calloc(strlen(prefix) + 10, 1);
		strcpy(str, prefix); strcat(str, ".bwt");  bwt_bwt[0] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".sa"); bwt_restore_sa(str, bwt_bwt[0]);
		strcpy(str, prefix); strcat(str, ".rbwt"); bwt_bwt[1] = bwt_restore_bwt(str);
		strcpy(str, prefix); strcat(str, ".rsa"); bwt_restore_sa(str, bwt_bwt[1]);
		free(str);
		bwt_bns = bns_restore(prefix);
	}
}

extern "C" void bwa_seed2genome_destroy()
{
	if (bwt_bwt[0])
		bwt_destroy(bwt_bwt[0]); 
	if (bwt_bwt[1])
		bwt_destroy(bwt_bwt[1]);
	if (bwt_bns)
		bns_destroy(bwt_bns) ;
}



extern "C" void bwa_seed2genome_numchr(int32_t *num)
{
	*num=bwt_bns->n_seqs;
	
}

extern "C" void bwa_seed2genome_descchr(int i, char*& desc)
{
	
	desc=bwt_bns->anns[i].name;
}

extern "C" void bwa_seed2genome_lenchr(int i, int32_t * len)
{
	*len=bwt_bns->anns[i].len;
	
}

extern "C" bwa_seq_t * bwa_seed2genome_map(const char* read, int read_len, int strand, uint64_t *num, uint64_t *sa_k, uint64_t *sa_l)
{
	bwa_seq_t *p=(bwa_seq_t*)calloc(1, sizeof(bwa_seq_t)) ;
	int n_seqs=1 ;
	
	int l = read_len ;
	//fprintf(stdout, "read=%s\n", read) ;
	
	p->tid = -1; // no assigned to a thread
	p->qual = NULL ;
	p->full_len = p->clip_len = p->len = l;
	p->seq = (ubyte_t*)calloc(p->len, 1);
	for (int i = 0; i != p->full_len; ++i)
	{
		p->seq[i] = nst_nt4_table[(int)read[i]];
		//fprintf(stdout, "seq[%i]=%i\n", i, p->seq[i]) ;
	}
	
	p->rseq = (ubyte_t*)calloc(p->full_len, 1);
	memcpy(p->rseq, p->seq, p->len);
	seq_reverse(p->len, p->seq, 0); // *IMPORTANT*: will be reversed back in bwa_refine_gapped()
	seq_reverse(p->len, p->rseq, bwt_opt->mode & BWA_MODE_COMPREAD);
	p->name = strdup("seq") ;
	p->cigar=NULL ;
	
	mybwa_cal_sa_reg_gap(0, bwt_bwt, n_seqs, p, bwt_opt);
	//fprintf(stdout, "n_aln=%i\n", p->n_aln) ;
	
	if (p->n_aln>0)
	{
		assert(p->n_aln<=2) ;
		if (p->aln[0].a==strand)
		{
			*sa_k=p->aln[0].k ;
			*sa_l=p->aln[0].l ;
			*num=*sa_l-*sa_k+1 ;
		} else
			if (p->n_aln>=2 && p->aln[1].a==strand)
			{
				*sa_k=p->aln[1].k ;
				*sa_l=p->aln[1].l ;
				*num=*sa_l-*sa_k+1 ;
			}
			else
			{
				*sa_k=1 ;
				*sa_l=0 ;
				*num=0 ;
			}
		//fprintf(stdout, "k=%lld, l=%lld\n", *sa_k, *sa_l) ;
	}
	else
	{
		*sa_k=1 ;
		*sa_l=0 ;
		*num=0 ;
	}

	/*p->sa = *sa_k ;
	p->c1 = 1 ;
	p->type=BWA_TYPE_UNIQUE ;
	
	mybwa_cal_pac_pos_core(bwt_bwt[0], bwt_bwt[1], p, 0, 0); 
				
	int len = pos_end(p) - p->pos; 
	int seq_id=0 ;
	
	bns_coor_pac2real(bwt_bns, p->pos, len, &seq_id) ;
	int pos = (int)(p->pos - bwt_bns->anns[seq_id].offset) ;
	fprintf(stdout, "seq_id=%i, pos=%i, n_aln=%i, multi=%i\n", seq_id, pos, p->n_aln, p->n_multi) ;
	*/
	
	return p ;
}

extern "C" void bwa_seed2genome_pos(uint64_t sa_pos, uint64_t *contig_id, uint64_t *contig_pos, bwa_seq_t *seq)
{
	bwa_seq_t *p=seq ;

	p->sa = sa_pos ;
	p->c1 = 1 ;
	p->type=BWA_TYPE_UNIQUE ;
	p->cigar=NULL ;
	p->strand=0 ;
	
	mybwa_cal_pac_pos_core(bwt_bwt[0], bwt_bwt[1], p, 0, 0); 
				
	uint64_t len = pos_end(p) - p->pos; 
	int seq_id=-1 ;
	
	bns_coor_pac2real(bwt_bns, p->pos, len, &seq_id) ;
	uint64_t pos = (int)(p->pos - bwt_bns->anns[seq_id].offset) ;
	
	if (false && sa_pos==461542)
	{
		fprintf(stdout, "seq_id=%i, pos=%lu, n_aln=%i, multi=%i, strand=%i\n", seq_id, pos, p->n_aln, p->n_multi, p->strand) ; 

		p->sa = 461542;//461970 ;
		p->c1 = 1 ;
		p->type=BWA_TYPE_UNIQUE ;
		p->cigar=NULL ;
		p->strand=1 ;
		
		mybwa_cal_pac_pos_core(bwt_bwt[0], bwt_bwt[1], p, 0, 0); 
		
		uint64_t len = pos_end(p) - p->pos; 
		int seq_id=-1 ;
		
		bns_coor_pac2real(bwt_bns, p->pos, len, &seq_id) ;
		uint64_t pos = (int)(p->pos - bwt_bns->anns[seq_id].offset) ;

		fprintf(stdout, "+++ seq_id=%i, pos=%lu, n_aln=%i, multi=%i, strand=%i\n", seq_id, pos, p->n_aln, p->n_multi, p->strand) ;
		//fprintf(stdout, "bwt->seq_len=%lld", (long long int)bwt_bwt[0]->seq_len) ;
		//fprintf(stdout, "reverse_bwt->seq_len=%lld", (long long int)bwt_bwt[1]->seq_len) ;
		
		//bwa_seq_t *a=NULL ;
		//fprintf(stdout, "error%lld", (long long int)a->sa) ;
	}

	*contig_id=seq_id ;
	*contig_pos=pos ;
}

extern "C" void bwa_seed2genome_cleanup_seq(bwa_seq_t *seq)
{
	bwa_free_read_seq(1, seq) ;
}

extern "C" void bwa_aln_my_core(const char *prefix, const char *read, gap_opt_t *opt)
{
	bwa_seed2genome_init(prefix, opt) ;

	uint64_t k=0, l=0, num=0 ;
	
	bwa_seq_t *seq = bwa_seed2genome_map(read, strlen(read), 0, &num, &k, &l) ;
	//fprintf(stdout, "k=%lu, l=%lu\n", k, l) ;
	
	for (int i=0; i<num; i++)
	{
		uint64_t contig_pos, contig_id ;
		bwa_seed2genome_pos(k+i, &contig_id, &contig_pos, seq) ;
	}
	bwa_seed2genome_cleanup_seq(seq) ;
	
	bwa_seed2genome_destroy() ;
	

	/*// initialization
	int i, n_seqs, tot_seqs = 0;
	bwa_seq_t *seqs;
	bwa_seqio_t *ks;
	ks = bwa_open_reads(bwt_opt->mode, fn_fa);

	// core loop
	//fwrite(opt, sizeof(gap_opt_t), 1, stdout);
	while ((seqs = bwa_read_seq(ks, 0x40000, &n_seqs, bwt_opt->mode & BWA_MODE_COMPREAD, bwt_opt->trim_qual)) != 0) {
		tot_seqs += n_seqs;

		mybwa_cal_sa_reg_gap(0, bwt_bwt, n_seqs, seqs, bwt_opt);

		for (i = 0; i < n_seqs; ++i) {
			bwa_seq_t *p = seqs + i;
			//bwa_aln2seq_core(p->n_aln, p->aln, p, 0, max_seed_num) ;
			//fprintf(stdout, "l=%i k=%i\n", p->aln->k, p->aln->l) ;
			
			for (uint64_t j=p->aln->k; j<= p->aln->l; j++)
			{
				p->sa = j ;
				p->c1=1 ;
				p->type=BWA_TYPE_UNIQUE ;
				
				mybwa_cal_pac_pos_core(bwt_bwt[0], bwt_bwt[1], &seqs[i], 0, 0); 
				
				int len = pos_end(p) - p->pos; // j is the length of the reference in the alignment
				int seq_id=0 ;
				
				bns_coor_pac2real(bwt_bns, p->pos, len, &seq_id) ;
				int pos = (int)(p->pos - bwt_bns->anns[seq_id].offset) ;
				
				fprintf(stdout, "seq_id=%i, pos=%i, n_aln=%i, multi=%i\n", seq_id, pos, p->n_aln, p->n_multi) ;
			}
		}

		bwa_free_read_seq(n_seqs, seqs);
	}

	// destroy
	bwa_seq_close(ks);
	*/
	
}

