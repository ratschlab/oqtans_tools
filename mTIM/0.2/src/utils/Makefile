include ../../bin/mTIM_config.sh

MEX=${MATLAB_MEX_PATH}
MKOCTFILE=${OCTAVE_MKOCT}
MATLAB_INCL=${MATLAB_INCLUDE_DIR}
SAMDIR=${SAMTOOLS_DIR}

all:	get_reads.mex interval_query.mex get_reads.mexa64 interval_query.mexa64
octave:	get_reads.mex interval_query.mex
matlab:	get_reads.mexa64 interval_query.mexa64

get_reads.mexa64:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	rm -f *.o
	${MEX} -g -O get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses -I$(MATLAB_INCL)

interval_query.mexa64:	get_string.c interval_query.c mmap.c find_interval.c
	rm -f *.o
	${MEX} -g -O interval_query.c get_string.c find_interval.c mmap.c -I$(MATLAB_INCL)

get_reads.mex:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	rm -f *.o
	${MKOCTFILE} -g --mex get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses

interval_query.mex: get_string.c interval_query.c mmap.c find_interval.c
	rm -f *.o
	${MKOCTFILE} -g --mex interval_query.c get_string.c find_interval.c mmap.c

clean:	
	rm -f *.o *.mexa64 *.mex 
