#!/bin/bash

. ../MasterTest.sh

CleanFiles wavelet.in wavelet.gnu cluster.gnu cluster.dat

TESTNAME='Wavelet analysis tests'
Requires netcdf

INPUT="-i wavelet.in"

cat > wavelet.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
rms @C,CA,N first

wavelet nb 10 s0 2 ds 0.25 type morlet correction 1.01 chival 0.25 \
        :1-22 out wavelet.gnu usemap
EOF
RunCpptraj "Wavelet analysis test"
DoTest wavelet.gnu.save wavelet.gnu

cat > wavelet.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
rms @C,CA,N first

wavelet nb 10 s0 2 ds 0.25 type morlet correction 1.01 chival 0.25 \
        :1-22 name DPDP \
        cluster clustermapout cluster.gnu clusterout cluster.dat \
          minpoints 66 epsilon 10.0
datafile cluster.gnu usemap palette kbvyw
EOF
RunCpptraj "Wavelet WAFEX test"
DoTest cluster.gnu.save cluster.gnu
DoTest cluster.dat.save cluster.dat

EndTest
exit 0
