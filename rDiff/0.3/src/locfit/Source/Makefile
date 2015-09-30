include ../../../bin/rdiff_config.sh

MEX=${RDIFF_MATLAB_MEX_PATH}
MKOCTFILE=${RDIFF_OCTAVE_MKOCT}
MATLAB_INCL=${RDIFF_MATLAB_INCLUDE_DIR}
LFSRC=liblfev.c liblocf.c libmut.c libtube.c
LIBDIR=../lib
INSTALLDIR=/usr/local
all: octave
octave:mexlf.mex mexpp.mex ../mex/mexlf.mex ../mex/mexpp.mex
matlab:mexlf mexpp

../mex/mexlf.mex: mexlf.mex
	cp mexlf.mex ../mex/
../mex/mexpp.mex: mexpp.mex
	cp mexpp.mex ../mex/

mexlf.mex: mexlf.c mlfut.c
	${MKOCTFILE} -v -g --mex $(MXFLAGS) $(CFLAGS) mexlf.c mlfut.c liblocf.c libmut.c libtube.c liblfev.c $(LDFLAGS)
mexpp.mex: mexpp.c mlfut.c
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
	rm -f *.o *.mex *.mexglx *.mexa64
FORCE:
