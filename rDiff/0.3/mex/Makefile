include ../bin/rdiff_config.sh

all:	get_reads.mex interval_overlap.mex get_reads.mexa64 interval_overlap.mexa64
octave:	get_reads.mex interval_overlap.mex
matlab:	get_reads.mexa64 interval_overlap.mexa64


get_reads.mexa64:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	rm -f *.o
	${MATLAB_MEX_PATH} -g -O get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${RDIFF_SAMTOOLS_INCLUDE_DIR} -L${RDIFF_SAMTOOLS_LIB_DIR} -lbam -lz -lcurses -I$(RDIFF_MATLAB_INCLUDE_DIR)

interval_overlap.mexa64:	interval_overlap.cpp
	${MATLAB_MEX_PATH} -g -O interval_overlap.cpp -I$(RDIFF_MATLAB_INCLUDE_DIR)

get_reads.mex:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	rm -f *.o
	${RDIFF_OCTAVE_MKOCT} -g --mex get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${RDIFF_SAMTOOLS_INCLUDE_DIR} -L${RDIFF_SAMTOOLS_LIB_DIR} -lbam -lz -lcurses

interval_overlap.mex: interval_overlap.cpp
	rm -f *.o
	${RDIFF_OCTAVE_MKOCT} -g --mex interval_overlap.cpp

clean:	
	rm -f *.o *.mexa64 *.mex 
