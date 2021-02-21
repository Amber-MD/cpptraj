include ../../config.h

OBJECTS=rngtest.o Random.o RNG.o RNG_Stdlib.o RNG_Marsaglia.o CpptrajStdio.o RNG_MersenneTwister.o RNG_PCG32.o xoshiro128plusplus.o RNG_Xoshiro128pp.o

all: rngtest

rngtest: $(OBJECTS)
	$(CXX) -o rngtest $(OBJECTS)

clean:
	/bin/rm -f $(OBJECTS) rngtest

rngtest.o: rngtest.cpp
	$(CXX) -c -o rngtest.o rngtest.cpp

Random.o: ../../src/Random.cpp
	$(CXX) -c -o Random.o ../../src/Random.cpp

RNG.o: ../../src/RNG.cpp
	$(CXX) -c -o RNG.o ../../src/RNG.cpp

RNG_Stdlib.o: ../../src/RNG_Stdlib.cpp
	$(CXX) -c -o RNG_Stdlib.o ../../src/RNG_Stdlib.cpp

RNG_Marsaglia.o: ../../src/RNG_Marsaglia.cpp
	$(CXX) -c -o RNG_Marsaglia.o ../../src/RNG_Marsaglia.cpp

RNG_MersenneTwister.o: ../../src/RNG_MersenneTwister.cpp
	$(CXX) $(DIRECTIVES) $(CXXFLAGS) -c -o RNG_MersenneTwister.o ../../src/RNG_MersenneTwister.cpp

RNG_PCG32.o: ../../src/RNG_PCG32.cpp
	$(CXX) -c -o RNG_PCG32.o ../../src/RNG_PCG32.cpp

CpptrajStdio.o: ../../src/CpptrajStdio.cpp
	$(CXX) -c -o CpptrajStdio.o ../../src/CpptrajStdio.cpp

RNG_Xoshiro128pp.o: ../../src/RNG_Xoshiro128pp.cpp
	$(CXX) -c -o RNG_Xoshiro128pp.o ../../src/RNG_Xoshiro128pp.cpp

xoshiro128plusplus.o: ../../src/xoshiro128plusplus.c
	$(CC) -c -o xoshiro128plusplus.o ../../src/xoshiro128plusplus.c
	
