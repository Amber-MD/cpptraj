#!/bin/bash

. ../MasterTest.sh

CleanFiles irspec.in irspec.dat raw.dat

TESTNAME='Infrared spectrum calculation test'
Requires netcdf maxthreads 10

INPUT='-i irspec.in'

cat > irspec.in <<EOF
parm ../Test_systemVF/systemVF.parm7
trajin ../Test_systemVF/systemVF.nc
infraredspec IR out irspec.dat maxlag 5 tstep 0.1 rawout raw.dat
EOF
RunCpptraj "$TESTNAME"
DoTest irspec.dat.save irspec.dat
DoTest raw.dat.save raw.dat

EndTest
exit 0
