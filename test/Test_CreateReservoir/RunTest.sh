#!/bin/bash

. ../MasterTest.sh

CleanFiles create.in testres.bins.nc
TESTNAME='Create RREMD NetCDF reservoir test'
Requires netcdf maxthreads 1

INPUT='-i create.in'

cat > create.in <<EOF
readdata Eamber.dat name ENE
readdata ../Test_ClusterDihedral/cvt.dat.save name CVT
parm ../tz2.parm7
trajin ../tz2.nc 1 10
createreservoir testres.bins.nc ene ENE:2 bin CVT:2 temp0 350 iseed 2112
EOF
RunCpptraj "$TESTNAME"
NcTest testres.bins.nc.save testres.bins.nc

EndTest
exit 0
