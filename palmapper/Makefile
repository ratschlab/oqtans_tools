SrcDir = ./src
OutDir = ./out
ObjDir = $(OutDir)/o
ExeDir = $(OutDir)/a

SVNVERSION = $(shell svnversion)

CC = g++
#CFLAGS = -Wall -ggdb -pg # debug 
#CFLAGS = -Wall -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow -O9 -fexpensive-optimizations -frerun-cse-after-loop -fcse-follow-jumps -finline-functions -fschedule-insns2 -fthread-jumps -fforce-addr -fstrength-reduce -funroll-loops -march=native -mtune=native -pthread # linux amd64 optimized
CFLAGS =  -O9 -Wall -ggdb -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow # generic
#CFLAGS = -O9 -Wall -g -pg -Wno-unused-parameter -Wformat -Wformat-security -Wimplicit -Wparentheses -Wshadow # generic
GMFLAGS = -DGM
INCLUDE =  -Ishogun/ -Idyn_prog/ -Isrc
LDFLAGS = 

SHOGUN_OBJ = $(ObjDir)/palmapper/shogun/init.o \
	$(ObjDir)/palmapper/shogun/Mathematics.o \
	$(ObjDir)/palmapper/shogun/io.o \
	$(ObjDir)/palmapper/shogun/Parallel.o \
	$(ObjDir)/palmapper/shogun/Version.o \
	$(ObjDir)/palmapper/shogun/SGObject.o \
	$(ObjDir)/palmapper/shogun/ShogunException.o \
	$(ObjDir)/palmapper/shogun/Signal.o

DYNPROG_OBJ = $(ObjDir)/palmapper/dyn_prog/Mathmatics_dp.o \
	$(ObjDir)/palmapper/dyn_prog/io_dp.o \
	$(ObjDir)/palmapper/dyn_prog/qpalma_dp.o \
	$(ObjDir)/palmapper/dyn_prog/debug_tools.o \
	$(ObjDir)/palmapper/dyn_prog/penalty_info_dp.o \
	$(ObjDir)/palmapper/dyn_prog/result_align.o \
	$(ObjDir)/palmapper/dyn_prog/fill_matrix.o

LANG_OBJ = $(ObjDir)/lang/Thread.o

PM_OBJ = $(ObjDir)/palmapper/GenomeMaps.o \
	$(ObjDir)/palmapper/QPalma.o \
	$(ObjDir)/palmapper/JunctionMap.o \
	$(ObjDir)/palmapper/align.o \
	$(ObjDir)/palmapper/TopAlignments.o \
	$(ObjDir)/palmapper/IntervalQuery.o \
	$(ObjDir)/palmapper/palmapper.o \
	$(ObjDir)/palmapper/print.o \
	$(ObjDir)/palmapper/Chromosome.o \
	$(ObjDir)/palmapper/Config.o \
	$(ObjDir)/palmapper/FileReporter.o \
	$(ObjDir)/palmapper/Genome.o \
	$(ObjDir)/palmapper/Hits.o \
	$(ObjDir)/palmapper/Mapper.o \
	$(ObjDir)/palmapper/QueryFile.o \
	$(ObjDir)/palmapper/Read.o \
	$(ObjDir)/palmapper/Statistics.o \
	$(ObjDir)/palmapper/Util.o \
	$(SHOGUN_OBJ) $(DYNPROG_OBJ) $(LANG_OBJ)

PMIDX_OBJ = $(ObjDir)/pmindex/init.o \
	$(ObjDir)/pmindex/printindex.o \
	$(ObjDir)/pmindex/usage.o \
	$(ObjDir)/pmindex/write.o \
	$(ObjDir)/pmindex/load.o \
	$(ObjDir)/pmindex/index.o \
	$(ObjDir)/pmindex/alloc.o \
	$(ObjDir)/pmindex/pmindex.o

CurrentDir := $(shell pwd)

all: palmapper pmindex

src/bwa/libbwa.a:
	@echo Compiling libbwa 
	(cd src/bwa && make libbwa.a)

src/bwa/bwa: src/bwa/libbwa.a
	@echo Compiling bwa
	(cd src/bwa && make bwa) 

bwa: src/bwa/bwa
	ln -sf src/bwa/bwa

palmapper: src/bwa/libbwa.a bwa $(PM_OBJ) src/palmapper/*.h 
	$(CC) $(CFLAGS) $(INCLUDE) $(PM_OBJ) $(LDFLAGS) -lpthread -lz -lm -Lsrc/bwa -lbwa -o palmapper
	ln -sf palmapper genomemapper

pmindex:  $(PMIDX_OBJ) src/pmindex/*.h src/pmindex/pmindex_symbols.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o pmindex $(PMIDX_OBJ)
	ln -sf pmindex gmindex 

clean:
	(cd src/bwa && make clean)
	rm -rf $(OutDir) palmapper pmindex genomemapper gmindex bwa

test:
	(cd testcase; make test)

release_pm:
	make clean; mkdir -p ../release; cd ..; rsync -av $(CurrentDir) release; cd release; rm -rf */.settings */.cproject */.project */.svn */*/.svn */*/*/.svn */*/*/*/.svn; tar czvf ../release.$(SVNVERSION).tar.gz .; cd ..; rm -rf release

#todo also delete shogun...
release_gm:
	make clean; mkdir -p ../release; cd ..; rsync -av $(CurrentDir) release; cd release; rm -rf */.settings */.cproject */.project */.svn */*/.svn */*/*/.svn; tar czvf ../genomemapper.$(SVNVERSION).tar.gz .; cd ..; rm -rf release

# generic rule for compiling c++
$(ObjDir)/%.o : $(SrcDir)/%.cpp $(SrcDir)/%.h
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/%.o : $(SrcDir)/%.cpp
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

$(ObjDir)/%.o : $(SrcDir)/%.c
	@echo Compiling $<
	@mkdir -p $(@D)
	@$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<
