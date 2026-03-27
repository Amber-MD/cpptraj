#!/bin/bash

. ../MasterTest.sh

CleanFiles min.in min.dat toface.dat test.toface.dat test.pdb

TESTNAME='Minimum non-self imaged distance test'
Requires netcdf maxthreads 10

INPUT="-i min.in"

cat > min.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
minimage m1 :1-13 :1-13 out min.dat
minimage m2 :1-13 :1-13 maskcenter out min.dat
minimage m3 :1-13 toface out toface.dat
EOF
RunCpptraj "Minimum Image test."
DoTest min.dat.save min.dat
DoTest toface.dat.save toface.dat

cat > test.pdb <<EOF
CRYST1    4.000    4.000    4.000  90.00  90.00  90.00 P 1           2         
MODEL     1 
ATOM      1  N   THR A   1S      2.000   2.000   2.000  1.00 65.41           N 
ENDMDL
MODEL     2 
ATOM      1  N   THR A   1S      1.000   2.000   2.000  1.00 65.41           N 
ENDMDL
MODEL     3 
ATOM      1  N   THR A   1S      3.000   2.000   2.000  1.00 65.41           N 
ENDMDL
MODEL     4 
ATOM      1  N   THR A   1S      2.000   1.000   2.000  1.00 65.41           N 
ENDMDL
MODEL     5 
ATOM      1  N   THR A   1S      2.000   3.000   2.000  1.00 65.41           N 
ENDMDL
MODEL     6 
ATOM      1  N   THR A   1S      2.000   2.000   1.000  1.00 65.41           N 
ENDMDL
MODEL     7 
ATOM      1  N   THR A   1S      2.000   2.000   3.000  1.00 65.41           N 
ENDMDL
END                                                                             
EOF
cat > min.in <<EOF
noprogress
parm test.pdb
trajin test.pdb
minimage Test toface @1 out test.toface.dat
run
EOF
RunCpptraj "Simple distance to unit cell face test."
DoTest test.toface.dat.save test.toface.dat

EndTest
exit 0
