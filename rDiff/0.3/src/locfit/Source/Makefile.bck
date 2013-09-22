include ../../../bin/rdiff_config.sh

MEX=${MATLAB_MEX_PATH}
MKOCTFILE=${OCTAVE_MKOCT}
MATLAB_INCL=${MATLAB_INCLUDE_DIR}
SAMDIR=${SAMTOOLS_DIR}
LFSRC=liblfev.c liblocf.c libmut.c libtube.c
LIBDIR=../lib
INSTALLDIR=/usr/local
#LDFLAGS=-L$(LIBDIR) #-Wl,-rpath $(INSTALLDIR)/lib 
#LDFLAGS= -lmut
all: mexlf mexpp
octave:mexlf_oct mexpp_oct
matlab:mexlf mexpp

mexlf_oct: mexlf.c mlfut.c
	${MKOCTFILE} -v -g --mex $(MXFLAGS) $(CFLAGS) mexlf.c mlfut.c liblocf.c libmut.c libtube.c liblfev.c $(LDFLAGS)
mexpp_oct: mexpp.c mlfut.c
	${MKOCTFILE} -v -g --mex $(MXFLAGS) $(CFLAGS) mexpp.c mlfut.c liblocf.c libmut.c libtube.c liblfev.c  $(LDFLAGS)
mexlf: mexlf.c mlfut.c
	${MEX} $(MXFLAGS) $(CFLAGS) mexlf.c mlfut.c $(LDFLAGS)
mexpp: mexpp.c mlfut.c
	${MEX} $(MXFLAGS) $(CFLAGS) mexpp.c mlfut.c $(LDFLAGS)
nodlls: mexlf.c mlfut.c mexpp.c $(LFSRC)
	${MEX} $(MXFLAGS) mexlf.c mlfut.c $(LFSRC)
	${MEX} $(MXFLAGS) mexpp.c mlfut.c $(LFSRC)
very-clean: clean
	rm -f mexlf.mexglx mexpp.mexglx
clean:
	rm -f *.o
FORCE:
	
