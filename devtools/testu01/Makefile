include ../../config.h

TESTU01_DIR=$(TESTU01_HOME)

#AR= ar cqs

#CXX=g++
CXXFLAGS += -I../../src

SOURCES=cpptraj_rng.cpp RNG_MersenneTwister.cpp RNG.cpp CpptrajStdio.cpp Random.cpp RNG_Stdlib.cpp RNG_Marsaglia.cpp RNG_PCG32.cpp RNG_Xoshiro128pp.cpp xoshiro128plusplus.cpp

OBJECTS=$(SOURCES:.cpp=.o)

all: a.out

main.o: main.c
	$(CC) -c -o main.o -Wall -O3 -I$(TESTU01_DIR)/include main.c

#cpptraj.a: $(OBJECTS)
#	-/bin/rm -f cpptraj.a
#	$(AR) cpptraj.a $(OBJECTS)

a.out: $(OBJECTS) main.o
	$(CXX) -o a.out main.o $(OBJECTS) -L$(TESTU01_DIR)/lib -ltestu01 -lprobdist -lmylib -lm
#	gcc -std=c99 -Wall -O3 -o a.out cpptraj.a -L$(TESTU01_DIR)/lib -ltestu01 -lprobdist -lmylib -lm

clean:
	-/bin/rm *.o a.out

