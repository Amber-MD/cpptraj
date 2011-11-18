include ../config.h

include cpptrajfiles

all: cpptraj$(SFX)

install: cpptraj$(SFX) 
	mv cpptraj$(SFX) $(CPPTRAJBIN) 

findDepend: FindDepend.o 
	$(CXX) -Wall -o findDepend FindDepend.o

depend: findDepend
	./findDepend $(SOURCE) > cpptrajdepend

dependclean:
	-/bin/rm FindDepend.o
	-/bin/rm findDepend

cpptraj$(SFX): $(OBJECTS) $(PTRAJ_OBJECTS) $(ARPACK) $(LAPACK) $(BLAS) 
	$(CXX) $(LDFLAGS) -o cpptraj$(SFX) $(OBJECTS) $(PTRAJ_OBJECTS) $(FLIBS_PTRAJ) $(NETCDFLIB) $(ZLIB) $(BZLIB) $(ARPACK) $(LAPACK) $(BLAS) $(FLIBS)


$(BLAS):
	cd ../../blas && $(MAKE) install

$(ARPACK):
	cd ../../arpack && $(MAKE) install

$(LAPACK):
	cd ../../lapack && $(MAKE) install 

thermo.o:  $(PTRAJ_HOME)/thermo.f
	$(FC) -c $(FFLAGS) -o $@ $(PTRAJ_HOME)/thermo.f

pubfft.o:  $(PTRAJ_HOME)/pubfft.f
	$(FC) -c $(FFLAGS) -o $@ $(PTRAJ_HOME)/pubfft.f

clean:
	/bin/rm -f $(OBJECTS) $(PTRAJ_OBJECTS) cpptraj$(SFX)

libclean:
	cd ../../blas && $(MAKE) clean
	cd ../../lapack && $(MAKE) clean
	cd ../../arpack && $(MAKE) clean
	/bin/rm -f $(BLAS) $(ARPACK) $(LAPACK)

# Header dependencies
include cpptrajdepend
