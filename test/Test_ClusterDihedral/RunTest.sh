#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in cd.dat cf.dat ci.dat 
RequiresNetcdf "clusterdihedral test"
RequiresNotParallel "clusterdihedral test"
# clusterdihedral
TOP="../tz2.parm7"
INPUT="ptraj.in"
cat > ptraj.in <<EOF
noprogress
trajin ../tz2.nc
clusterdihedral out cd.dat framefile cf.dat clusterinfo ci.dat :2-12 phibins 2 psibins 2
EOF
RunCpptraj "clusterdihedral test."
DoTest cd.dat.save cd.dat
DoTest cf.dat.save cf.dat
DoTest ci.dat.save ci.dat

EndTest

exit 0
