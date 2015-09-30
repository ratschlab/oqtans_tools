// Authors: Korbinian Schneeberger and Joerg Hagmann
// Copyright (C) 2008 by Max-Planck Institute for Developmental Biology, Tuebingen, Germany

#include "Hits.h"
#include "align.h"
#include "Genome.h"

double WORST_SCORE;
double WORST_MM_SCORE;

int check_mm(Read const &read, Chromosome const &chr, int genome_pos, int readpos, int ori, char conversion)
{
	if (_config.BSSEQ) {
		char diff = (ori == 1? chr[genome_pos] : get_compl_base(chr[genome_pos])) - read.data()[readpos];

		if ( (conversion == 1 && diff != 0 && diff != -17) ||
		     (conversion == 2 && diff != 0 && diff !=   6) ||
		     !unique_base(read.data()[readpos]) )
		{
			return 1;
		} else {
			return 0;
		}
	}
	else {
		if ( (ori == 1  && (chr[genome_pos] != read.data()[readpos])) ||
		(ori == -1 && get_compl_base(chr[genome_pos]) != read.data()[readpos]) ||
		!(unique_base(read.data()[readpos]))
		)
		{
			return 1;
		} else {
			return 0;
		}
	}
}

// returns 1 if alignment is gapless, 0 if there are gaps and -1 if nr of allowed MMs is exceeded
int kbound_overhang_alignment(Read const &read, HIT* hit, int offset, int readstart, int start, int end, unsigned short int hitreadpos, Chromosome const &chromosome, char orientation, unsigned char mismatches)
{
	// global vars -> local vars
	char const * const READ = read.data();
	int Read_length = read.length();
	int K = _config.NUM_GAPS;
	int Max_read_length = _config.MAX_READ_LENGTH;
	double Mismatch_score = _config.MM_SCORE;
	double Match_score = _config.M_SCORE;
	double Gap_score = _config.GAP_SCORE;
	double Worst_score = WORST_SCORE;
	double Worst_mm_score = WORST_MM_SCORE;

	//if (_config.STATISTICS) _stats.NUM_ALIGNMENTS++;

	int length;
	int chrstart;
	char offset_comp = 0;

	EDIT_OPS edit_op[Config::MAX_EDIT_OPS];
	//memcpy(edit_op, hit->edit_op, mismatches * sizeof(EDIT_OPS));
	for (int ii=0; ii<mismatches; ii++)
		edit_op[ii] = hit->edit_op[ii] ;

	int i, j, h;
	// TODO: dd make this stack vars?
	double** M = (double**) malloc ((2 * K + 1) * sizeof(double*));
	char** T = (char**) malloc ((2 * K + 1) * sizeof(char*));

	if (M == NULL || T == NULL) {
		fprintf(stderr, "[kbound_overhang_alignment] ERROR Could not allocate memory\n");
		exit(1);
	}

	for (i=0; i!=2*K+1; ++i) {
		*(M+i) = (double*) malloc ((Max_read_length + K * 2) * sizeof(double));
		*(T+i) = (char*) malloc ((Max_read_length + K * 2) * sizeof(char));
		if (*(M+i) == NULL || *(T+i) == NULL) {
			fprintf(stderr, "[kbound_overhang_alignment] ERROR Could not allocate memory\n");
			exit(1);
		}
	}

	if (readstart == 0) {
		length = hitreadpos + offset - 1;
		offset_comp = offset;
		if (orientation == '+') chrstart = start - length - 1;
			else chrstart = end + length - 1;
	}
	else {
		length = Read_length - (end - start) - hitreadpos + offset;
		if (orientation == '+') chrstart = end;
			else chrstart = start - 2;
	}

	// Initialization:
	double score_before = mismatches * Mismatch_score;
	if (readstart == 0) {
		for (i=0; i!=2*K+1; ++i) {
			M[i][offset] = score_before;
			//T[i][offset] = '0';
		}
		j = (K-offset < 0)? 0: K-offset;
		for (i=0; i!=j; ++i) {
			M[0][i+1+offset] = score_before + (i+1) * Gap_score;
			//T[0][i+1+offset] = UP;
		}
	}
	else {
		for (i = 0; i <= K; ++i) {
			M[i][0] = score_before + i * Gap_score;
			if (i!=0) {
				M[0][i] = score_before + i * Gap_score;

				//T[i][0] = LEFT;
				//T[0][i] = UP;
			}
		}
		//T[0][0] = '0';
	}

	// Alignment:
	int c, min_i = -1;
	char best_score_on_diag = 0;
	double score;
	double best_score = Worst_score + 1;
	double column_score = best_score;

	for (i = 1; i <= length; ++i) {

		for (h = -K; h <= K; ++h) {
			j = i + h;

				if (j <= length && j >= 1) {
					if (j>K) c = j - K;
						else c = 0;

				if ( !(readstart == 0 && j <= offset) && !(readstart != 0 && j > length-offset) ) {

					//if (_config.STATISTICS) ++_stats.CELLS_OVERHANG;

					// Score-Function:
					if (orientation == '+') {
						if (chrstart+i-1>=(int)chromosome.length())
						{
							for (i=0; i!=2*K+1; ++i) {
								free(M[i]);
								free(T[i]);
							}
							free(M);
							free(T);
							return -1;
						}

						if (_config.BSSEQ) {
							char diff = chromosome[chrstart+i-1] - READ[readstart+j-offset_comp-1];
							if ( (hit->conversion == 1 && diff != 0 && diff != -17) ||
							     (hit->conversion == 2 && diff != 0 && diff !=   6) ||
							     !unique_base(READ[readstart+j-offset_comp-1]) )
								score = Mismatch_score;
							else	score = Match_score;
						} else {
							if (chromosome[chrstart+i-1] != READ[readstart+j-offset_comp-1] || !unique_base(READ[readstart+j-offset_comp-1]))
								score = Mismatch_score;
							else	score = Match_score;
						}
					}
					else {
						if (chrstart-i+1<0)
						{
							for (i=0; i!=2*K+1; ++i) {
								free(M[i]);
								free(T[i]);
							}
							free(M);
							free(T);
							return -1;
						}

						if (_config.BSSEQ) {
							char diff = get_compl_base(chromosome[chrstart+i-1]) - READ[readstart+j-offset_comp-1];
							if ( (hit->conversion == 1 && diff != 0 && diff != -17) ||
							     (hit->conversion == 2 && diff != 0 && diff !=   6) ||
							     !unique_base(READ[readstart+j-offset_comp-1]) )
								score = Mismatch_score;
							else	score = Match_score;
						} else {
							if (get_compl_base(chromosome[chrstart-i+1]) != READ[readstart+j-offset_comp-1] || !unique_base(READ[readstart+j-offset_comp-1]))
								score = Mismatch_score;
							else	score = Match_score;
						}
					}

					M[i-c][j] = M[i-c-(j<=K)][j-1] + score;
					//T[i-c][j] = DIAGONAL;	// traceback diagonally

					if ((i-j-1 <= K) && (i-j-1 >= -K)) {
						if (M[i-c-1][j] + Gap_score <= M[i-c][j]) {

							if (i == j) {
								//if (_config.STATISTICS) _stats.GAPS_ENCOUNTERED[0]++;
								// cleanup:
								for (i=0; i!=2*K+1; ++i) {
									free(M[i]);
									free(T[i]);
								}
								free(M);
								free(T);
								return 0;
							}

							M[i-c][j] = M[i-c-1][j] + Gap_score;
							//T[i-c][j] = LEFT;	// traceback to the left

						}
					}
					if ((i-j+1 <= K) && (i-j+1 >= -K)) {
						if (M[i-c+(j>K)][j-1] + Gap_score <= M[i-c][j]) {

							if (i == j) {
								//if (_config.STATISTICS) _stats.GAPS_ENCOUNTERED[0]++;
								// cleanup:
								for (i=0; i!=2*K+1; ++i) {
									free(M[i]);
									free(T[i]);
								}
								free(M);
								free(T);
								return 0;
							}

							M[i-c][j] = M[i-c+(j>K)][j-1] + Gap_score;
							//T[i-c][j] = UP;	// traceback to upper

						}
					}

					if (readstart != 0 && j == length - offset) {
						if ( (!best_score_on_diag && M[i-c][j] < best_score) || (best_score_on_diag && M[i-c][j] <= best_score) ) {
							// ensures that if score of most right bottom cell is equal to another cell in the same row, this other cell is reported with the best score
							best_score = M[i-c][j];
							min_i = i;
							if (i == j) best_score_on_diag = 1;
						}
					}

					if (M[i-c][j] < column_score) column_score = M[i-c][j];

				}
			}

		} //for h

		if (column_score > Worst_score) {

			// readstart==0: score in most right bottom cell cannot become less than worst score
			// best_score > WORST_SCORE: there is no cell in the bottom row which has a valid alignment score
			if (readstart == 0 || best_score > Worst_score) {
				//if (_config.STATISTICS) _stats.TOO_MANY_MMS[0]++;
				// cleanup:
				for (i=0; i!=2*K+1; ++i) {
					free(M[i]);
					free(T[i]);
				}
				free(M);
				free(T);
				return -1;
			}
			// if best_score has been found before most right bottom cell was processed
			/*else if (i < length - offset) {	// readstart here is != 0
				_stats.GAPS_ENCOUNTERED[1]++;
				return 0;
			}*/

			break;

		} else {
			column_score = Worst_score + 1;
		}

	} //for i


	if (readstart == 0) {
		j = length;
		i = K - (j<K) * (K - j);
		best_score = M[i][j];
		min_i = j;
	}
	else {
		j = length-offset;
		i = min_i - (j>K) * (j - K);
	}

	if (best_score > Worst_score) {
		//if (_config.STATISTICS) _stats.TOO_MANY_MMS[1]++;
		// cleanup:
		for (i=0; i!=2*K+1; ++i) {
			free(M[i]);
			free(T[i]);
		}
		free(M);
		free(T);
		return -1;
	}

	// if best score is worse than max nr of allowed mismatches...
	if (best_score > Worst_mm_score) {
		// ...and best score is in the most right bottom cell, there is no better path with gaps -> so we have too many MMs
		if (min_i == j) {
			//if (_config.STATISTICS) _stats.TOO_MANY_MMS[1]++;
			// cleanup:
			for (i=0; i!=2*K+1; ++i) {
				free(M[i]);
				free(T[i]);
			}
			free(M);
			free(T);
			return -1;
		}
		// ..and there is another cell with this score other than the most right bottom cell -> there is a path with a valid score, BUT with gaps!
		else {
			//if (_config.STATISTICS) _stats.GAPS_ENCOUNTERED[2]++;
			// cleanup:
			for (i=0; i!=2*K+1; ++i) {
				free(M[i]);
				free(T[i]);
			}
			free(M);
			free(T);
			return 0;
		}
	}

	// FOR READSTART != 0 ONLY: if best path doesn't start from the most right bottom corner of the matrix -> gaps are in best alignment, return 0
	//if (min_i != j && M[i][j] != M[K-(j<K)*(K-j)][j]) {	@TODO!!!!!!!
	if (min_i != j) {	// even if score in most right bottom cell is equal to best score -> path with gaps needs less edit ops -> perform global algnm.
		//if (_config.STATISTICS) _stats.GAPS_ENCOUNTERED[2]++;
		// cleanup:
		for (i=0; i!=2*K+1; ++i) {
			free(M[i]);
			free(T[i]);
		}
		free(M);
		free(T);
		return 0;
	}
	// FOR READSTART == 0: traceback has to start in most right bottom corner



	// Traceback (if there had been gaps the procedure returned 0, so only mismatches are left -> only diagonal traceback)
	int readpos;
	if (readstart == 0) readpos = j - offset;
		else readpos = readstart + j;
	i = K - (j<K) * (K - j);

	int mms = (int) (M[i][j] / Mismatch_score) - (int) (score_before / Mismatch_score);


	while (j != (readstart==0)*offset && mms != 0) {
		if (M[i][j] != M[i-(j<=K)][j-1]) {
			// create mismatch:

			if (orientation == '+') edit_op[mismatches].pos = readpos;
				else		edit_op[mismatches].pos = Read_length - readpos + 1;
			edit_op[mismatches].mm = 1;
			mismatches++;
			//assert(mismatches<Config::MAX_EDIT_OPS) ;

			--mms;
		}
		--readpos;
		--j;
		if (j<K) --i;
	}


	// cleanup:
	for (i=0; i!=2*K+1; ++i) {
		free(M[i]);
		free(T[i]);
	}
	free(M);
	free(T);

	// write in hit-structure:
	//memcpy(hit->edit_op, edit_op, mismatches * sizeof(EDIT_OPS));
	for (int ii=0; ii<mismatches; ii++)
		hit->edit_op[ii] = edit_op[ii] ;
	hit->mismatches = mismatches;
	assert(mismatches<Config::MAX_EDIT_OPS) ;

	// successfully aligned
	return 1;
}


// k-bound global alignment algorithm:
int kbound_global_alignment(Read const &read, HIT* hit, unsigned short int hitreadpos, unsigned int start, unsigned int end, Chromosome const &chromosome, char orientation, int Num_edit_ops)
{
	// global vars -> local vars
	char const * const READ = read.data();
	int Read_length = read.length();
	unsigned int Chr_length = chromosome.length();
	double Gap_score = _config.GAP_SCORE;
	double Mismatch_score = _config.MM_SCORE;
	double Match_score = _config.M_SCORE;
	double Worst_score = WORST_SCORE;
	char Gaps_most_right = _config.GAPS_MOST_RIGHT;
	int Num_mismatches = _config.NUM_MISMATCHES;
	char Stringent_gaplimit = _config.STRINGENT_GAPLIMIT;
	int K = _config.NUM_GAPS;
	int Max_read_length = _config.MAX_READ_LENGTH;

	// arguments:
	unsigned char mismatches = 0;
	unsigned char gaps = 0;
	signed char start_offset = 0;
	signed char end_offset = 0;

	int i, j, h;

	// initialize alignment and traceback matrices and mismatch array
	double** M = (double**) malloc ((2 * K + 1) * sizeof(double*));
	char** T = (char**) malloc ((2 * K + 1) * sizeof(char*));

	if (M == NULL || T == NULL) {
		fprintf(stderr, "[kbound_global_alignment] ERROR Could not allocate memory\n");
		exit(1);
	}

	for (i=0; i!=2*K+1; ++i) {
		*(M+i) = (double*) malloc ((Max_read_length + K * 2) * sizeof(double));
		*(T+i) = (char*) malloc ((Max_read_length + K * 2) * sizeof(char));
		if (*(M+i) == NULL || *(T+i) == NULL) {
			fprintf(stderr, "[kbound_global_alignment] ERROR Could not allocate memory\n");
			exit(1);
		}
	}

	EDIT_OPS edit_op[Config::MAX_EDIT_OPS];

	int chrstart, chrend ;
	int offset_front, offset_end;
	if (orientation == '+') {
		chrstart = start - hitreadpos;			// 0-initialized
		if ((int)chrstart < K) 
			offset_front = chrstart;
		else offset_front = K;

		chrend = chrstart + Read_length;		// 1-initialized
		if (chrend + K > (int)Chr_length) 
			offset_end = Chr_length - chrend;
		else offset_end = K;

		chrstart -= offset_front;
	}
	else {
		chrstart = end + hitreadpos - 2;		// 0-initialized
		if (chrstart + K >= (int)Chr_length) offset_front = Chr_length - chrstart - 1;
			else offset_front = K;

		chrend = chrstart - Read_length + 1;		// 0-initialized
		if (chrend - K < 0) offset_end = chrend;
			else offset_end = K;

		chrstart += offset_front;			// 0-initialized
	}

	int length = Read_length + offset_front + offset_end;



	// Initialization:
	for (i=0; i!=2*K+1-(K-offset_front); ++i) {

		M[i][offset_front] = 0;
		T[i][offset_front] = '0';

		j = (K-offset_front < 0)? 0: K-offset_front;
		for (h=0; h!=j; ++h) {
			M[0][h+1+offset_front] = (h+1) * Gap_score;
			T[0][h+1+offset_front] = UP;
		}
	}


	// Alignment:
	int c;
	double best_score = Worst_score + 1;
	double column_score = best_score;
	double score;
	int min_i = 0;	// start of traceback
	char best_score_on_diag = 0;

	for (i = 1; i <= length; ++i) {

		for (h = -K; h <= K; ++h) {
			j = i + h;

				if (j <= length && j >= 1) {
					if (j>K) c = j - K;
						else c = 0;

				if ( j > offset_front && j <= length-offset_end ) 
				{

//					if (_config.STATISTICS) ++_stats.CELLS_GLOBAL;

					// Score-Function:
					if (orientation == '+') {
						if (chrstart+i-1>=(int)chromosome.length())
						{
							for (i=0; i!=2*K+1; ++i) {
								free(M[i]);
								free(T[i]);
							}
							free(M);
							free(T);
							return -1;
						}

						if (_config.BSSEQ) {
							int diff = chromosome[chrstart+i-1] - READ[j-offset_front-1];
							if ( (hit->conversion == 1 && diff != 0 && diff != -17) ||
							     (hit->conversion == 2 && diff != 0 && diff !=   6) ||
							     !unique_base(READ[j-offset_front-1]) )
								score = _config.MM_SCORE;
							else	score = _config.M_SCORE;
							
						} else {
							if (chromosome[chrstart+i-1] != READ[j-offset_front-1] || !unique_base(READ[j-offset_front-1]))
								score = _config.MM_SCORE;
							else	score = _config.M_SCORE;
						}
					}
					else {
						if (chrstart-i+1<0)
						{
							for (i=0; i!=2*K+1; ++i) {
								free(M[i]);
								free(T[i]);
							}
							free(M);
							free(T);
							return -1;
						}
						
						if (_config.BSSEQ) {
							int diff = get_compl_base(chromosome[chrstart-i+1]) - READ[j-offset_front-1];
							if ( (hit->conversion == 1 && diff != 0 && diff != -17) ||
							     (hit->conversion == 2 && diff != 0 && diff !=   6) ||
							     !unique_base(READ[j-offset_front-1]) )
								score = _config.MM_SCORE;
							else	score = _config.M_SCORE;
							
						} else {
							if (get_compl_base(chromosome[chrstart-i+1]) != READ[j-offset_front-1] || !unique_base(READ[j-offset_front-1]))
								score = Mismatch_score;
							else	score = Match_score;
						}
					}

					M[i-c][j] = M[i-c-(j<=K)][j-1] + score;
					T[i-c][j] = 'D';	// traceback diagonally

					if ((i-j+1 <= K) && (i-j+1 >= -K)) {
						// gap in chr
						if (M[i-c+(j>K)][j-1] + Gap_score == M[i-c][j] &&
								((orientation == '+' && Gaps_most_right) || (orientation == '-' && !Gaps_most_right))) {
							M[i-c][j] = M[i-c+(j>K)][j-1] + Gap_score;
							T[i-c][j] = 'U';	// traceback to upper with gaps most right
						}
						else if (M[i-c+(j>K)][j-1] + Gap_score < M[i-c][j]) {
							M[i-c][j] = M[i-c+(j>K)][j-1] + Gap_score;
							T[i-c][j] = 'U';	// traceback to upper
						}
					}

					if ((i-j-1 <= K) && (i-j-1 >= -K)) {
						// gap in read
						if (M[i-c-1][j] + Gap_score == M[i-c][j] &&
								((orientation == '+' && Gaps_most_right) || (orientation == '-' && !Gaps_most_right))) {
							M[i-c][j] = M[i-c-1][j] + Gap_score;
							T[i-c][j] = 'L';	// traceback to the left with gaps most right
						}
						else if (M[i-c-1][j] + Gap_score < M[i-c][j]) {
							M[i-c][j] = M[i-c-1][j] + Gap_score;
							T[i-c][j] = 'L';	// traceback to the left
						}
					}

					// Remember best score, i.e. start of traceback
					if (j == length - offset_end) {

						// gaps in reads preferred to gaps in chr:
						if ( (i <= j && M[i-c][j] <= best_score) || (i > j && ((best_score_on_diag && M[i-c][j] <  best_score)
							|| (!best_score_on_diag && M[i-c][j] <= best_score))) ) {

							best_score = M[i-c][j];
							min_i = i;
							if (i == j) best_score_on_diag = 1;
						}
					}
					// best_score preference:
					//	1. diagonal
					//	2. i > j (gaps in reads)
					//	3. i < j (gaps in chr)

					if (M[i-c][j] < column_score) column_score = M[i-c][j];

				}
			}

		} //for h

		if (column_score > Worst_score) {
			if (best_score > Worst_score) {
				//if (_config.STATISTICS) _stats.BREAK_GLOBAL_ALIGNMENT[0]++;
				// cleanup:
				for (i=0; i!=2*K+1; ++i) {
					free(M[i]);
					free(T[i]);
				}
				free(M);
				free(T);
				return -1;
			}
			else {
				//if (_config.STATISTICS) _stats.BREAK_GLOBAL_ALIGNMENT[1]++;
				break;
			}
		} else {
			column_score = Worst_score + 1;
		}

	} //for i

	if (best_score > Worst_score) {
                //if (_config.STATISTICS) _stats.BREAK_GLOBAL_ALIGNMENT[0]++;
        // cleanup:
		for (i=0; i!=2*K+1; ++i) {
			free(M[i]);
			free(T[i]);
		}
		free(M);
		free(T);
		return -1;
	}


	j = length - offset_end;
	i = min_i - (j>K)*(j-K);
	int chrpos;
	if (orientation == '+') {
		chrpos = chrstart + min_i - 1;
		end_offset += min_i - j;	// correcting the read pos on chr for later output
	}
	else {
		chrpos = chrstart - min_i + 1;
		start_offset -= min_i - j;	// correcting the read pos on chr for later output
	}
	int readpos = Read_length - 1;




	// Traceback:

	while (T[i][j] != '0' && M[i][j] != 0) {

		switch (T[i][j]) {
			case 'D': {

				//if (CHR_SEQ[hit->chromosome][chrpos] != READ[readpos] || !unique_base(READ[readpos])) {
				//if (M[i-(j<=K)][j-1] == M[i][j] - _config.MM_SCORE) {					// doesn't work for doubles!!
				//if (fabs(M[i-(j<=K)][j-1] - M[i][j] + _config.MM_SCORE) > 1E-9) {		// works, but requires math.h!
				if (M[i-(j<=K)][j-1] != M[i][j]) {
					if ((mismatches-gaps) < Num_mismatches && mismatches < Num_edit_ops) {
						if (orientation == '+') edit_op[mismatches].pos = readpos + 1;
							else		edit_op[mismatches].pos = Read_length - readpos;
						edit_op[mismatches].mm = 1;
						mismatches++;
						//assert(mismatches<Config::MAX_EDIT_OPS) ;
					}
					else {
						//if (_config.STATISTICS) _stats.BREAK_TB_IN_GLOBAL_ALIGNMENT++;
						// cleanup:
						for (i=0; i!=2*K+1; ++i) {
							free(M[i]);
							free(T[i]);
						}
						free(M);
						free(T);
						return -1;	// discard hit
					}
				}

				i = i-(j<=K);
				j--;
				if (orientation == '+') chrpos--;
					else chrpos++;
				readpos--;
				break;
			}

			case 'L': {
				if (mismatches < Num_edit_ops && (!Stringent_gaplimit || gaps < K)) {
					//if (_config.STATISTICS && hit->gaps >= _config.NUM_GAPS) W++;
					if (orientation == '+') edit_op[mismatches].pos = readpos + 2;
						else		edit_op[mismatches].pos = Read_length - readpos;
					edit_op[mismatches].mm = 0;
					mismatches++;
					//assert(mismatches<Config::MAX_EDIT_OPS) ;
					gaps++;
				}
				else {
					//if (_config.STATISTICS) _stats.BREAK_TB_IN_GLOBAL_ALIGNMENT++;
					// cleanup:
					for (i=0; i!=2*K+1; ++i) {
						free(M[i]);
						free(T[i]);
					}
					free(M);
					free(T);
					return -1;	// discard hit
				}

				i--;
				if (orientation == '+') chrpos--;
					else chrpos++;
				break;
			}

			case 'U': {

				if (mismatches < Num_edit_ops && (!Stringent_gaplimit || gaps < K)) {
					//if (_config.STATISTICS && hit->gaps >= _config.NUM_GAPS) W++;
					if (orientation == '+') edit_op[mismatches].pos = -readpos - 1;
						else		edit_op[mismatches].pos = -Read_length + readpos;
					edit_op[mismatches].mm = 0;
					mismatches++;
					//assert(mismatches<Config::MAX_EDIT_OPS) ;
					gaps++;
				}
				else {
					//if (_config.STATISTICS) _stats.BREAK_TB_IN_GLOBAL_ALIGNMENT++;
					// cleanup:
					for (i=0; i!=2*K+1; ++i) {
						free(M[i]);
						free(T[i]);
					}
					free(M);
					free(T);
					return -1;	// discard hit
				}

				i = i+(j>K);
				j--;
				readpos--;
				break;
			}
		}

	}

	if (orientation == '+') start_offset += i + (j>K)*(j-K) - j;
		else end_offset -= i + (j>K)*(j-K) - j;

	// write in hit-structure:
	hit->mismatches = mismatches;
	assert(mismatches<=Config::MAX_EDIT_OPS) ;
	hit->gaps = gaps;
	hit->start_offset = start_offset;
	hit->end_offset = end_offset;
	//memcpy(hit->edit_op, edit_op, mismatches * sizeof(EDIT_OPS));
	for (int ii=0; ii<mismatches; ii++)
		hit->edit_op[ii]=edit_op[ii] ;

	// cleanup:
	for (i=0; i!=2*K+1; ++i) {
		free(M[i]);
		free(T[i]);
	}
	free(M);
	free(T);

	// successfully aligned
	return mismatches;
}
