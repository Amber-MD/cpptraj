#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles general.in distance.dat rmsd.dat rmsda.dat phi2.dat PhiPsi.dat \
           test.crd a1.dat Restart/* Restart test.nc r4.dat a2.dat.gz \
           a3.dat.bz2 r2.dat r3-nofit.dat

TESTNAME='General tests'
ERR=0
NotParallel "$TESTNAME"
((ERR = ERR + $?))
# Check libraries
CheckNetcdf "General tests"
((ERR = ERR + $?))

if [ $ERR -ne 0 ] ; then
  SkipTest "$TESTNAME"
  exit 0
fi

CheckZlib
CheckBzlib

if [ ! -e 'Restart' ] ; then
  mkdir Restart
fi

cat > general.in <<EOF
noprogress
distance d1 :1 :2 out distance.dat
distance d2 :1 :22 out distance.dat
trajin ../tz2.pdb
trajin ../tz2.crd.gz
trajout test.crd
trajout test.nc netcdf
trajout Restart/test.rst7 restart
reference ../tz2.rst7
angle a1 :2@CA :3@CA :4@CA out a1.dat
angle a2 :2@CA :3@CA :4@CA out a2.dat.gz
angle a3 :2@CA :3@CA :4@CA out a3.dat.bz2
#rmsd r4 first out r4.dat :2-10 nofit
rmsd r1 ref tz2.rst7 out rmsd.dat
rmsd r1a refindex 0 out rmsda.dat
rmsd r2 :1-4@C,CA,N ref tz2.rst7 :1-4@C,CA,N out r2.dat
rmsd r3 :2-5@CA ref tz2.rst7 :2-5@CA out r3-nofit.dat nofit
dihedral dh1 :1@C :2@N :2@CA :2@C out phi2.dat
#dihedral dh1 :2@C :3@N :3@CA :3@C
dihedral phi3 :2@C :3@N :3@CA :3@C out PhiPsi.dat
dihedral psi3 :3@N :3@CA :3@C :4@N out PhiPsi.dat
datafile PhiPsi.dat noxcol 
parm ../tz2.parm7
parm ../tz2.ortho.parm7
parm ../DPDP.parm7
reference ../DPDP.nc parm DPDP.parm7
trajin ../tz2.ortho.nc parm tz2.ortho.parm7
rms r4 ref DPDP.nc :2-11 nofit out r4.dat
#  trajin ../tz2.nc
    # Skip this
trajin ../DPDP.nc parm ../DPDP.parm7
trajin ../tz2.crd   \
 2 \
3 \
2

EOF

INPUT="general.in"
TOP="../tz2.parm7"
RunCpptraj "General tests"

DoTest distance.dat.save distance.dat
DoTest rmsd.dat.save rmsd.dat
DoTest rmsda.dat.save rmsda.dat
DoTest phi2.dat.save phi2.dat
DoTest PhiPsi.dat.save PhiPsi.dat
DoTest test.crd.save test.crd
DoTest a1.dat.save a1.dat
DoTest test.rst7.213.save Restart/test.rst7.213
NcTest test.nc.save test.nc
DoTest r4.dat.save r4.dat
# NOTE: a2.dat.gz comparison allowed to fail on windows; differences caused
#       by different newline characters in compressed file. Macs also seem to
#       occasionally fail this test, even though decompressed diffs are the same
if [ "$CPPTRAJ_TEST_OS" = 'linux' ] ; then
  DoTest a2.dat.gz.save a2.dat.gz
fi
DoTest a3.dat.bz2.save a3.dat.bz2
DoTest r2.dat.save r2.dat
DoTest r3-nofit.dat.save r3-nofit.dat

EndTest

exit 0
