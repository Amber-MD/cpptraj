include ../config.h

include cpptrajfiles

LIBDIR=/opt/Amber/GIT/amber/lib
FLIBS_PTRAJ=../../ptraj/thermo.o ../../ptraj/pubfft.o $(LIBDIR)/arpack.a $(LIBDIR)/lapack.a $(LIBDIR)/blas.a   -L/cottus/opt/intel/Compiler/11.1/072/lib/intel64/ -lifport -lifcore -lsvml

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

cpptraj$(SFX): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o cpptraj$(SFX) $(OBJECTS) $(FLIBS_PTRAJ) $(NETCDFLIB) $(ZLIB) $(BZLIB)

clean:
	/bin/rm -f $(OBJECTS) cpptraj$(SFX)

# Header dependencies
include cpptrajdepend
