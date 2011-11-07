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

cpptraj$(SFX): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o cpptraj$(SFX) $(OBJECTS) $(NETCDFLIB) $(ZLIB) $(BZLIB)

clean:
	/bin/rm -f $(OBJECTS) cpptraj$(SFX)

# Test for evec.o
evec.o: ../../ptraj/evec.o
	/bin/cp ../../ptraj/evec.o .

# Header dependencies
include cpptrajdepend
