SRCS= Mathmatics_dp.cpp\
		fill_matrix.cpp\
		qpalma_dp.cpp\
		result_align.cpp\
		debug_tools.cpp\
		penalty_info_dp.cpp\
		print_align.cpp\
		io_dp.cpp

HDRS= Mathmatics_dp.h\
		common_dp.h\
		config_dp.h\
		debug_tools.h\
		fill_matrix.h\
		io_dp.h\
		penalty_info_dp.h\
		qpalma_dp.h

OBJS = $(SRCS:%.cpp=%.o)

#CXXFLAGS=-O3 -fPIC
#CXXFLAGS=-O3 -fPIC -pg -fprofile-arcs

#CXXFLAGS=-O3 -fPIC -ggdb -pg
#CXXFLAGS=-O3 -fPIC -ggdb -lprofiler #-fno-inline
CXXFLAGS= -O -fPIC -g -ggdb 

LDFLAGS=-lprofiler -L/fml/ag-raetsch/home/fabio/own_libs/libunwind/lib -lunwind-x86_64 -lunwind

all: $(OBJS) $(HDRS)
	g++ $(CXXFLAGS) -shared $(OBJS) -o libqpalma.so

test: $(OBJS) $(HDRS) 
	g++ $(CXXFLAGS) $(LDFLAGS) -o test_fm  debug_tools.o Mathmatics.o io.o penalty_info.o fill_matrix.o test_fill_matrix.cpp 

clean:
	rm *.o *.so

