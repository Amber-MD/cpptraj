#!/bin/bash

. ../MasterTest.sh

if [[ ! -e Restart ]] ; then
  mkdir Restart
fi

# Clean
CleanFiles distance.dat rmsd.dat rmsda.dat phi2.dat PhiPsi.dat test.crd a1.dat Restart/* test.nc r4.dat a2.dat.gz a3.dat.bz2 r2.dat r3-nofit.dat

INPUT="general.in"
TOP="../trpcage.parm7"
RunCpptraj "General tests"

DoTest distance.dat.save distance.dat
DoTest rmsd.dat.save rmsd.dat
DoTest rmsda.dat.save rmsda.dat
DoTest phi2.dat.save phi2.dat
DoTest PhiPsi.dat.save PhiPsi.dat
DoTest test.crd.save test.crd
DoTest a1.dat.save a1.dat
DoTest test.rst7.102.save Restart/test.rst7.102
DoTest test.nc.save test.nc
DoTest r4.dat.save r4.dat
DoTest a2.dat.gz.save a2.dat.gz
DoTest a3.dat.bz2.save a3.dat.bz2
DoTest r2.dat.save r2.dat
DoTest r3-nofit.dat.save r3-nofit.dat
CheckTest

EndTest

exit 0
