include ../../config.h

OBJECTS=rngtest.o Random.o RNG.o RNG_Stdlib.o RNG_Marsaglia.o CpptrajStdio.o RNG_MersenneTwister.o

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

CpptrajStdio.o: ../../src/CpptrajStdio.cpp
	$(CXX) -c -o CpptrajStdio.o ../../src/CpptrajStdio.cpp
