include ../config.h

include cpptrajfiles

OBJECTS=$(SOURCES:.cpp=.o) $(CSOURCES:.c=.o) 

all: cpptraj$(SFX)

install: cpptraj$(SFX) 
	mv cpptraj$(SFX) $(CPPTRAJBIN) 

findDepend: FindDepend.o 
	$(CXX) -Wall -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCES) $(CSOURCES) > cpptrajdepend

dependclean:
	-/bin/rm FindDepend.o
	-/bin/rm findDepend

cpptraj$(SFX): $(OBJECTS) pubfft.o $(ARPACK) $(LAPACK) $(BLAS) 
	$(CXX) $(LDFLAGS) -o cpptraj$(SFX) $(OBJECTS) pubfft.o $(NETCDFLIB) $(ZLIB) $(BZLIB) $(ARPACK) $(LAPACK) $(BLAS) $(FLIBS)


$(BLAS): $(BLAS_HOME)/Makefile
	cd $(BLAS_HOME) && $(MAKE) install

$(ARPACK): $(ARPACK_HOME)/Makefile
	cd $(ARPACK_HOME) && $(MAKE) install

$(LAPACK): $(LAPACK_HOME)/Makefile
	cd $(LAPACK_HOME) && $(MAKE) install

pubfft.o:  pubfft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pubfft.F90

clean:
	/bin/rm -f $(OBJECTS) pubfft.o cpptraj$(SFX)

libclean: $(BLAS_HOME)/Makefile $(ARPACK_HOME)/Makefile $(LAPACK_HOME)/Makefile
	cd $(BLAS_HOME) && $(MAKE) clean
	cd $(LAPACK_HOME) && $(MAKE) clean
	cd $(ARPACK_HOME) && $(MAKE) clean
	/bin/rm -f $(BLAS) $(ARPACK) $(LAPACK)

# Header dependencies
include cpptrajdepend
