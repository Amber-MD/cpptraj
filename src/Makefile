# Cpptraj standalone Makefile
include ../config.h

include cpptrajfiles

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o) $(TIMER)

include ambpdbfiles

AMBPDB_OBJECTS=$(AMBPDBSOURCES:.cpp=.o)

all: cpptraj$(SFX)

install: cpptraj$(SFX) ambpdb$(SFX) nproc$(SFX) 
	mv cpptraj$(SFX) $(CPPTRAJBIN)
	(if [ -e "ambpdb$(SFX)" ]; then \
	   mv ambpdb$(SFX) $(CPPTRAJBIN); \
	 fi;) 

findDepend: FindDepend.o 
	$(CXX) -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) $(CSOURCES) AmbPDB.cpp > cpptrajdepend

dependclean:
	/bin/rm -f FindDepend.o findDepend

libcpptraj: $(LIBCPPTRAJ)

libcpptraj$(SHARED_SUFFIX): $(OBJECTS) $(FFT_DEPEND) $(READLINE_TARGET)
	$(CXX)  -shared -o libcpptraj$(SHARED_SUFFIX) $(OBJECTS) $(FFT_LIB) $(CPPTRAJ_LIB) $(LDFLAGS) $(READLINE)
	/bin/mv libcpptraj$(SHARED_SUFFIX) $(CPPTRAJLIB)

nolibcpptraj:
	@(echo "Error: Cannot build libcpptraj; re-configure with '-shared'" ; exit 1 ; )

cpptraj$(SFX): $(OBJECTS) $(FFT_DEPEND) $(READLINE_TARGET)
	$(CXX)  -o cpptraj$(SFX) $(OBJECTS) $(READLINE) $(FFT_LIB) $(CPPTRAJ_LIB) $(LDFLAGS)

ambpdb: $(AMBPDB_OBJECTS)
	$(CXX) -o ambpdb $(AMBPDB_OBJECTS) $(LDFLAGS)

ambpdb.exe: $(AMBPDB_OBJECTS)
	$(CXX) -o ambpdb.exe $(AMBPDB_OBJECTS) $(LDFLAGS)

ambpdb.OMP:
	@echo ""
	@echo "INFO: ambpdb is not OpenMP-enabled. Skipping."
	@echo ""

ambpdb.MPI:
	@echo ""
	@echo "INFO: ambpdb is not MPI-enabled. Skipping."
	@echo ""

nproc: ;

nproc.OMP: ;

nproc.MPI:
	cd ../util && $(MAKE)

nproc.exe: ;

readline/libreadline.a: $(READLINE_HOME)/Makefile
	cd $(READLINE_HOME) && $(MAKE) all

noreadline:
	@echo ""
	@echo "Skipping READLINE support"
	@echo ""

$(FFT_DEPEND):  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

ReadLine.o: ReadLine.cpp
	$(CXX) $(CXXFLAGS) -I$(READLINE_HOME) -c -o $@ ReadLine.cpp

clean:
	/bin/rm -f $(OBJECTS) $(FFT_DEPEND) cpptraj$(SFX) AmbPDB.o ambpdb
	cd $(READLINE_HOME) && $(MAKE) clean

uninstall:
	/bin/rm -f $(CPPTRAJBIN)/cpptraj $(CPPTRAJBIN)/cpp_ambpdb readline/libreadline.a $(CPPTRAJLIB)/libcpptraj$(SHARED_SUFFIX)

# Header dependencies
include cpptrajdepend
