#!/bin/bash

. ../MasterTest.sh

TESTNAME='Random number generator tests'

CleanFiles rng.in random.dat mt.dat

INPUT='-i rng.in'

UNITNAME='Test Marsaglia, Stdlib, PCG32, and Xoshiro128++ RNGs'
cat > rng.in <<EOF
rng setdefault marsaglia createset Marsaglia    settype int count 10 seed 10 out random.dat
rng setdefault stdlib    createset Stdlib       settype int count 10 seed 10 out random.dat
rng setdefault pcg32     createset PCG32        settype int count 10 seed 10 out random.dat
rng setdefault xo128     createset Xoshiro128++ settype int count 10 seed 10 out random.dat
list
EOF
RunCpptraj "$UNITNAME"

UNITNAME='Test Mersenne Twister RNG'
CheckFor c++11
if [ $? -eq 0 ] ; then
  cat > rng.in <<EOF
rng setdefault mt createset MT    settype int count 10 seed 10 out mt.dat
EOF
  RunCpptraj "$UNITNAME"
fi

EndTest
