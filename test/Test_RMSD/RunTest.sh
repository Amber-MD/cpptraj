#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat rms.mass.in rmsd.mass.dat rms.reftraj.in \
           rmsd.reftraj.dat tz2.norotate.crd tz2.rotate.crd rmatrices.dat \
           rmsd.refcoords.dat rms.dat NoMod.dat NoMod.crd.save NoMod.crd \
           vecs.dat Previous.dat

TESTNAME='RMSD tests'
Requires netcdf

# Test rmsd, mass-weighted rmsd, rmsd to reference traj.
UNITNAME='Basic RMSD tests'
CheckFor maxthreads 10
if [ $? -eq 0 ] ; then
  TOP='../tz2.truncoct.parm7'
  INPUT='rms.in'
  cat > rms.in <<EOF
noprogress
trajin ../tz2.truncoct.nc

rms Res2-11 first :2-11 out rmsd.dat
rms Res2-11_mass first :2-11 out rmsd.mass.dat mass \
  savevectors combined vecsout vecs.dat
rms Res2_11_traj reftraj ../tz2.truncoct.nc :2-11 out rmsd.reftraj.dat
run
removedata Res2_11_traj
loadtraj ../tz2.truncoct.nc name TZ2
rms Res2_11_traj reftraj TZ2 :2-11 out rmsd.refcoords.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest rmsd.dat.save rmsd.dat
  DoTest rmsd.mass.dat.save rmsd.mass.dat
  DoTest vecs.dat.save vecs.dat -a 0.000001
  DoTest rmsd.reftraj.dat.save rmsd.reftraj.dat
  DoTest rmsd.reftraj.dat.save rmsd.refcoords.dat
fi

# Test RMS rotate/norotate, generation of rotation matrices
UNITNAME='RMS coordinate rotation/rotation matrices test'
CheckFor maxthreads 10
if [ $? -eq 0 ] ; then
  TOP=''
  INPUT='-i rms.in'
  cat > rms.in <<EOF
parm ../tz2.parm7 [NOWAT] 
reference ../tz2.nc parm [NOWAT] 1 [first] 
parm ../tz2.truncoct.parm7 [WAT]
trajin ../tz2.truncoct.nc parm [WAT]
strip :WAT
rms NOROT ref [first] norotate @CA
outtraj tz2.norotate.crd parm [WAT]
rms ROT ref [first] out rms.dat @CA savematrices matricesout rmatrices.dat
outtraj tz2.rotate.crd parm [WAT]
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.norotate.crd.save tz2.norotate.crd
  DoTest tz2.rotate.crd.save tz2.rotate.crd
  DoTest rmatrices.dat.save rmatrices.dat
fi

# Test RMS nomod
TOP=''
INPUT="-i rms.in"
cat > rms.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
outtraj NoMod.crd.save
rms First_CA :2-12@CA out NoMod.dat nomod
trajout NoMod.crd
EOF
RunCpptraj "RMS fit with no coordinates modification test."
DoTest NoMod.dat.save NoMod.dat
DoTest NoMod.crd.save NoMod.crd

# Test RMS 'previous'
UNITNAME='RMS fit to previous test'
CheckFor notparallel
if [ $? -eq 0 ] ; then
  TOP=''
  INPUT='-i rms.in'
  cat > rms.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
rms ToPrevious :2-12@CA previous out Previous.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest Previous.dat.save Previous.dat
fi

EndTest

exit 0
