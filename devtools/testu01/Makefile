TESTU01_DIR=/home/droe/Programs/TestU01/1.2.3

SOURCES=main.cpp ../../src/RNG_MersenneTwister.cpp ../../src/RNG.cpp ../../src/CpptrajStdio.cpp ../../src/Random.cpp ../../src/RNG_Stdlib.cpp ../../src/RNG_Marsaglia.cpp ../../src/RNG_PCG32.cpp ../../src/RNG_Xoshiro128pp.cpp ../../src/xoshiro128plusplus.cpp

all: a.out

a.out: $(SOURCES) 
	g++ -o a.out -I$(TESTU01_DIR)/include -L$(TESTU01_DIR)/lib $(SOURCES) -ltestu01 -lprobdist -lmylib

clean:
	-/bin/rm *.o a.out

