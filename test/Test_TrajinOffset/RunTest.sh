#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles cpptraj.offset.in rem.crd.combined gzip.rem.crd.combined \
           bzip2.rem.crd.combined phi.dat.save phi.dat

TESTNAME='Trajectory read with offset tests'
Requires maxthreads 13
# Test 1
cat > cpptraj.offset.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7
trajin rem.crd.000 1 10 2 
trajin rem.crd.001 2 10 3
trajin rem.crd.002 3 8 2
trajin rem.crd.003 1 10 5
trajout rem.crd.combined 
EOF
INPUT="-i cpptraj.offset.in"
RunCpptraj "Normal trajectory read with offsets."
DoTest rem.crd.save rem.crd.combined 

# Test 2
UNITNAME='Gzipped trajectory read with offsets'
CheckFor zlib
if [ $? -eq 0 ] ; then
  cat > cpptraj.offset.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7
trajin rem.crd.000.gz 1 10 2 
trajin rem.crd.001.gz 2 10 3
trajin rem.crd.002.gz 3 8 2
trajin rem.crd.003.gz 1 10 5
trajout gzip.rem.crd.combined 
EOF
  INPUT="-i cpptraj.offset.in"
  RunCpptraj "$UNITNAME"
  DoTest rem.crd.save gzip.rem.crd.combined
fi

# Test 3
UNITNAME='Bzip2ed trajectory read with offsets'
CheckFor bzlib
if [ $? -eq 0 ] ; then
  cat > cpptraj.offset.in <<EOF
noprogress
parm ala2.99sb.mbondi2.parm7
trajin rem.crd.000.bz2 1 10 2 
trajin rem.crd.001.bz2 2 10 3
trajin rem.crd.002.bz2 3 8 2
trajin rem.crd.003.bz2 1 10 5
trajout bzip2.rem.crd.combined 
EOF
  INPUT="-i cpptraj.offset.in"
  RunCpptraj "$UNITNAME"
  DoTest rem.crd.save bzip2.rem.crd.combined
fi

# Test 4
UNITNAME='Start argument offset'
CheckFor maxthreads 6
if [ $? -eq 0 ] ; then
  cat > cpptraj.offset.in <<EOF
parm ala2.99sb.mbondi2.parm7
trajin rem.crd.000 5 10
dihedral phi0 :1@C :2@N :2@CA :2@C out phi.dat.save noheader
run
clear trajin
trajin rem.crd.000 -5
dihedral phi1 :1@C :2@N :2@CA :2@C out phi.dat noheader
EOF
  RunCpptraj "$UNITNAME"
  DoTest phi.dat.save phi.dat
fi

EndTest

exit 0
