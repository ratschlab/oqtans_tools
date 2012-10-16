CXX = g++
CXXFLAGS = -O9 -g -Wall -fmessage-length=0 -ULARGEFILE64_SOURCE
CFLAGS= 

MATLAB_DIR = /fml/ag-raetsch/share/software

MEX = $(MATLAB_DIR)/matlabR2011a/bin/mex ${CFLAGS}

MATLAB_INCL = -I$(MATLAB_DIR)/matlabR2011a/extern/include

SAMDIR = /fml/ag-raetsch/share/software/samtools

get_reads.mexa64:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	rm -f *.o
	${MEX} -g get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses $(MATLAB_INCL)

get_reads.mex:	get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp
	# mkoctfile --mex get_reads.cpp get_reads_direct.cpp mex_input.cpp read.cpp -I${SAMDIR} -L${SAMDIR} -lbam -lz -lcurses


interval_query.mexa64: get_string.c interval_query.c mmap.c find_interval.c
	rm -f *.o
	${MEX} -g -O interval_query.c get_string.c find_interval.c mmap.c -I/fml/ag-raetsch/share/software/matlab/extern/include

interval_query.mex: get_string.c interval_query.c mmap.c find_interval.c
	rm -f *.o
	# mkoctfile -g --mex interval_query.c get_string.c find_interval.c mmap.c


all:	get_reads.mexa64 get_reads.mex interval_query.mexa64 interval_query.mex

clean:	
	rm -f *.o *.mexa64 *.mex get_reads_direct
