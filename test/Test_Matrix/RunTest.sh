#!/bin/bash

. ../MasterTest.sh

CleanFiles matrix.in mtest.dat.save mtest.*.dat evecs.10.dat \
           tz2.dist.ca.matrix.dat.save tz2.dist.ca.matrix.dat

TESTNAME='Matrix Tests'
Requires maxthreads 10

INPUT="-i matrix.in"
cat > matrix.in <<EOF
parm 1rrb_vac.prmtop
trajin 1rrb_vac.mdcrd
matrix dist @CA out mtest.0.dat byres
matrix dist @N @CA out mtest.1.dat byres
matrix dist @CA out mtest.2.dat bymask
matrix dist @N @CA out mtest.3.dat bymask
matrix correl @N @C out mtest.4.dat
matrix covar @N @C out mtest.5.dat
matrix mwcovar @N @C out mtest.6.dat
matrix dist @CA out mtest.7.dat
matrix idea @CA out mtest.8.dat
matrix correl @CA out mtest.9.dat
matrix covar @CA out mtest.10.dat
matrix mwcovar @CA out mtest.11.dat
matrix dist @N @C out mtest.12.dat 
matrix distcovar :1-4@CA out mtest.13.dat
EOF
RunCpptraj "$TESTNAME"
DoTest mtest.dat.0.save mtest.0.dat
DoTest mtest.dat.1.save mtest.1.dat
DoTest mtest.dat.2.save mtest.2.dat
DoTest mtest.dat.3.save mtest.3.dat
DoTest mtest.dat.4.save mtest.4.dat
DoTest mtest.dat.5.save mtest.5.dat
DoTest mtest.dat.6.save mtest.6.dat
DoTest mtest.dat.7.save mtest.7.dat
DoTest mtest.dat.8.save mtest.8.dat
DoTest mtest.dat.9.save mtest.9.dat
DoTest mtest.dat.10.save mtest.10.dat
DoTest mtest.dat.11.save mtest.11.dat
DoTest mtest.dat.12.save mtest.12.dat
DoTest mtest.dat.13.save mtest.13.dat

# Test reading symmetric matrix
# NOTE: Currently disabled due to eigenvector sign flips causing false test errors.
ReadSymmMatrix() {
cat > matrix.in <<EOF
readdata mtest.dat.10.save name MyMatrix read2d
diagmatrix MyMatrix out evecs.10.dat vecs 3
EOF
RunCpptraj "Read symmetric matrix data test."
DoTest evecs.10.dat.save evecs.10.dat
}

# Test start/stop/offset args
cat > matrix.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 5 100 10
matrix dist @CA out tz2.dist.ca.matrix.dat.save
EOF
RunCpptraj "Generate matrix with trajectory start/stop/offset"
cat > matrix.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
matrix dist @CA out tz2.dist.ca.matrix.dat start 5 stop 100 offset 10
EOF
RunCpptraj "Generate matrix with action stop/start/offset"
DoTest tz2.dist.ca.matrix.dat.save tz2.dist.ca.matrix.dat

EndTest
  
exit 0
