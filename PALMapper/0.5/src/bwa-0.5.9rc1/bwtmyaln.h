
extern "C" void bwa_seed2genome_init(const char *prefix, gap_opt_t *opt);
extern "C" void bwa_seed2genome_destroy() ;
extern "C" bwa_seq_t * bwa_seed2genome_map(const char* read, int read_len, int strand, uint64_t *num, uint64_t *sa_k, uint64_t *sa_l) ;
extern "C" void bwa_seed2genome_pos(uint64_t sa_pos, uint64_t *contig_id, uint64_t *contig_pos, bwa_seq_t *seq) ;
extern "C" void bwa_seed2genome_cleanup_seq(bwa_seq_t *seq) ;
extern "C" void bwa_seed2genome_numchr(int32_t *num);
extern "C" void bwa_seed2genome_descchr(int i, char *&desc);
extern "C" void bwa_seed2genome_lenchr(int i, int32_t *len);








