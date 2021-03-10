#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in random.crd around.overlap.rst7 around.rst7 overlap.rst7 norestrict.rst7

TESTNAME='randomizeions test'
Requires maxthreads 1 zlib

INPUT="ptraj.in"
TOP="adh206.ff10.tip3p.parm7.gz"
cat > ptraj.in <<EOF
trajin adh206.tip3p.rst7.gz
rng setdefault marsaglia
randomizeions @Na+ around :1-16 by 5.0 overlap 3.0 seed 113698 originalalgorithm
trajout random.crd title "Test"
EOF
RunCpptraj "randomizeions test (original algorithm)"
DoTest random.crd.save random.crd

# Around and overlap
cat > ptraj.in <<EOF
trajin adh206.tip3p.rst7.gz
rng setdefault marsaglia
randomizeions :Na+ around :1-16 by 5.0 overlap 3.0 seed 1
trajout around.overlap.rst7
EOF
RunCpptraj "Randomize (around, overlap)"
DoTest around.overlap.rst7.save around.overlap.rst7

# Around 
cat > ptraj.in <<EOF
trajin adh206.tip3p.rst7.gz
rng setdefault marsaglia
randomizeions :Na+ around :1-16 by 5.0 allowoverlap seed 1
trajout around.rst7
EOF
RunCpptraj "Randomize (around)"
DoTest around.rst7.save around.rst7

# Overlap
cat > ptraj.in <<EOF
trajin adh206.tip3p.rst7.gz
rng setdefault marsaglia
randomizeions :Na+ overlap 3.0 seed 1
trajout overlap.rst7
EOF
RunCpptraj "Randomize (overlap)"
DoTest overlap.rst7.save overlap.rst7

# No restrictions
cat > ptraj.in <<EOF
trajin adh206.tip3p.rst7.gz
rng setdefault marsaglia
randomizeions :Na+ allowoverlap seed 1
trajout norestrict.rst7
EOF
RunCpptraj "Randomize (no restrictions)"
DoTest norestrict.rst7.save norestrict.rst7

EndTest

exit 0
