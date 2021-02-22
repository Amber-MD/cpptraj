#!/bin/bash

. ../MasterTest.sh

TESTNAME='Random number generator tests'

CleanFiles rng.in random.dat

INPUT='-i rng.in'

UNITNAME='Test Marsaglia RNG'
cat > rng.in <<EOF
rng setdefault marsaglia createset Marsaglia    settype int count 10 seed 10 out random.dat
rng setdefault stdlib    createset Stdlib       settype int count 10 seed 10 out random.dat
rng setdefault pcg32     createset PCG32        settype int count 10 seed 10 out random.dat
rng setdefault xo128     createset Xoshiro128++ settype int count 10 seed 10 out random.dat
list
EOF
RunCpptraj "$UNITNAME"

EndTest
