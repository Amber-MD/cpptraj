include ../../config.h

#OBJECTS=rngtest.o Random.o RNG.o RNG_Stdlib.o RNG_Marsaglia.o CpptrajStdio.o RNG_MersenneTwister.o RNG_PCG32.o xoshiro128plusplus.o RNG_Xoshiro128pp.o
SOURCES=rngtest.cpp ../../src/RNG_MersenneTwister.cpp ../../src/RNG.cpp ../../src/CpptrajStdio.cpp ../../src/Random.cpp ../../src/RNG_Stdlib.cpp ../../src/RNG_Marsaglia.cpp ../../src/RNG_PCG32.cpp ../../src/RNG_Xoshiro128pp.cpp ../../src/xoshiro128plusplus.cpp

all: rngtest

debug:
	make CXXFLAGS='-O0 -g' CFLAGS='-O0 -g' rngtest

rngtest: $(SOURCES) 
	$(CXX) $(DIRECTIVES) $(CXXFLAGS) -o rngtest $(SOURCES) 

clean:
	/bin/rm -f  rngtest

rngtest.o: rngtest.cpp
	$(CXX) $(CXXFLAGS) -c -o rngtest.o rngtest.cpp

Random.o: ../../src/Random.cpp
	$(CXX) $(CXXFLAGS) -c -o Random.o ../../src/Random.cpp

RNG.o: ../../src/RNG.cpp
	$(CXX) $(CXXFLAGS) -c -o RNG.o ../../src/RNG.cpp

RNG_Stdlib.o: ../../src/RNG_Stdlib.cpp
	$(CXX) $(CXXFLAGS) -c -o RNG_Stdlib.o ../../src/RNG_Stdlib.cpp

RNG_Marsaglia.o: ../../src/RNG_Marsaglia.cpp
	$(CXX) $(CXXFLAGS) -c -o RNG_Marsaglia.o ../../src/RNG_Marsaglia.cpp

RNG_MersenneTwister.o: ../../src/RNG_MersenneTwister.cpp
	$(CXX) $(DIRECTIVES) $(CXXFLAGS) -c -o RNG_MersenneTwister.o ../../src/RNG_MersenneTwister.cpp

RNG_PCG32.o: ../../src/RNG_PCG32.cpp
	$(CXX) $(CXXFLAGS) -c -o RNG_PCG32.o ../../src/RNG_PCG32.cpp

CpptrajStdio.o: ../../src/CpptrajStdio.cpp
	$(CXX) $(CXXFLAGS) -c -o CpptrajStdio.o ../../src/CpptrajStdio.cpp

RNG_Xoshiro128pp.o: ../../src/RNG_Xoshiro128pp.cpp
	$(CXX) $(CXXFLAGS) -c -o RNG_Xoshiro128pp.o ../../src/RNG_Xoshiro128pp.cpp

xoshiro128plusplus.o: ../../src/xoshiro128plusplus.cpp
	$(CXX) $(CXXFLAGS) -c -o xoshiro128plusplus.o ../../src/xoshiro128plusplus.cpp
	
