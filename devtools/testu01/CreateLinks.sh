#!/bin/bash

for FILE in ../../src/RNG_MersenneTwister.cpp ../../src/RNG.cpp ../../src/CpptrajStdio.cpp ../../src/Random.cpp ../../src/RNG_Stdlib.cpp ../../src/RNG_Marsaglia.cpp ../../src/RNG_PCG32.cpp ../../src/RNG_Xoshiro128pp.cpp ../../src/xoshiro128plusplus.cpp ; do
  ln -s $FILE .
done

