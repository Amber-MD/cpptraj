#!/bin/bash

. ../MasterTest.sh

TESTNAME='Random number generator tests'

CleanFiles rng.in random.dat mt.dat stdlib.dat pcg32.dat

INPUT='-i rng.in'

UNITNAME='Test Marsaglia, PCG32, and Xoshiro128++ RNGs'
cat > rng.in <<EOF
rng setdefault marsaglia createset Marsaglia    settype int count 10 seed 10 out random.dat
rng setdefault xo128     createset Xoshiro128++ settype int count 10 seed 10 out random.dat
list
EOF
RunCpptraj "$UNITNAME"
DoTest random.dat.save random.dat

UNITNAME='PCG32 RNG test.'
# PCG32 does not yet compile for windows.
CheckFor notos windows
if [ $? -eq 0 ] ; then
  cat > rng.in <<EOF
rng setdefault pcg32     createset PCG32        settype int count 10 seed 10 out pcg32.dat
list
EOF
  RunCpptraj "$UNITNAME"
  DoTest pcg32.dat.save pcg32.dat
fi

UNITNAME='Test Stdlib RNG'
# Stdlib rand implementation is different on windows/osx
CheckFor testos Linux
if [ $? -eq 0 ] ; then
  cat > rng.in <<EOF
rng setdefault stdlib    createset Stdlib       settype int count 10 seed 10 out stdlib.dat
list
EOF
  RunCpptraj "$UNITNAME"
  DoTest stdlib.dat.save stdlib.dat
fi

UNITNAME='Test Mersenne Twister RNG'
CheckFor c++11
if [ $? -eq 0 ] ; then
  cat > rng.in <<EOF
rng setdefault mt createset MT    settype int count 10 seed 10 out mt.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest mt.dat.save mt.dat
fi

EndTest
