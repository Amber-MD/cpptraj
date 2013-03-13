# Cpptraj stand-alone make file
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
	/bin/rm -f FindDepend.o findDepend

cpptraj$(SFX): $(OBJECTS) pub_fft.o $(READLINE)
	$(CXX) -o cpptraj$(SFX) $(OBJECTS) pub_fft.o $(READLINE) $(LDFLAGS)


$(READLINE): $(READLINE_HOME)/Makefile
	cd $(READLINE_HOME) && $(MAKE)

pub_fft.o:  pub_fft.F90
	$(FC) $(FPPFLAGS) -c $(FFLAGS) -o $@ pub_fft.F90

ReadLine.o: ReadLine.cpp
	$(CXX) $(CXXFLAGS) -I$(READLINE_HOME) -c -o $@ ReadLine.cpp

clean:
	/bin/rm -f $(OBJECTS) pub_fft.o cpptraj$(SFX)
	cd $(READLINE_HOME) && $(MAKE) clean

# Header dependencies
include cpptrajdepend
