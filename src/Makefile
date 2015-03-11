# Cpptraj standalone Makefile
include ../config.h

include cpptrajfiles

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o) $(TIMER)

include ambpdbfiles

AMBPDB_OBJECTS=$(AMBPDBSOURCES:.cpp=.o) MpiRoutines.o

all: cpptraj$(SFX)

install: cpptraj$(SFX) ambpdb$(SFX) 
	mv cpptraj$(SFX) $(CPPTRAJBIN)
	(if [ -e "ambpdb" ]; then \
	   mv ambpdb $(CPPTRAJBIN); \
	 fi;) 

findDepend: FindDepend.o 
	$(CXX) -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) $(CSOURCES) AmbPDB.cpp > cpptrajdepend

dependclean:
	/bin/rm -f FindDepend.o findDepend

libcpptraj: checkshared $(OBJECTS) pub_fft.o $(READLINE)
	$(CXX)  -shared -o libcpptraj.so $(OBJECTS) pub_fft.o $(LDFLAGS) $(READLINE)
	/bin/mv libcpptraj.so $(CPPTRAJLIB)

cpptraj$(SFX): $(OBJECTS) pub_fft.o $(READLINE)
	$(CXX)  -o cpptraj$(SFX) $(OBJECTS) pub_fft.o $(READLINE) $(LDFLAGS)


ambpdb: $(AMBPDB_OBJECTS)
	$(CXX) -o ambpdb $(AMBPDB_OBJECTS) $(LDFLAGS)

ambpdb.OMP:
	@echo ""
	@echo "INFO: ambpdb is not OpenMP-enabled. Skipping."
	@echo ""

ambpdb.MPI:
	@echo ""
	@echo "INFO: ambpdb is not MPI-enabled. Skipping."
	@echo ""

$(READLINE): $(READLINE_HOME)/Makefile
	cd $(READLINE_HOME) && $(MAKE) all

pub_fft.o:  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

ReadLine.o: ReadLine.cpp
	$(CXX) $(CXXFLAGS) -I$(READLINE_HOME) -c -o $@ ReadLine.cpp

clean:
	/bin/rm -f $(OBJECTS) pub_fft.o cpptraj$(SFX) AmbPDB.o ambpdb
	cd $(READLINE_HOME) && $(MAKE) clean

uninstall:
	/bin/rm -f $(CPPTRAJBIN)/cpptraj $(CPPTRAJBIN)/cpp_ambpdb readline/libreadline.a $(CPPTRAJLIB)/libcpptraj.so

# Header dependencies
include cpptrajdepend
