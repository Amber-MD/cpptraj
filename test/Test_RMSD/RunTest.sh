#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd.dat rms.mass.in rmsd.mass.dat rms.reftraj.in \
           rmsd.reftraj.dat tz2.norotate.crd tz2.rotate.crd rmatrices.dat \
           rmsd.refcoords.dat

CheckNetcdf
TOP="../tz2.truncoct.parm7"
INPUT="rms.in"

# Test rmsd, mass-weighted rmsd, rmsd to reference traj.
cat > rms.in <<EOF
noprogress
trajin ../tz2.truncoct.nc

rms Res2-11 first :2-11 out rmsd.dat
rms Res2-11_mass first :2-11 out rmsd.mass.dat mass
rms Res2_11_traj reftraj ../tz2.truncoct.nc :2-11 out rmsd.reftraj.dat
run
removedata Res2_11_traj
loadtraj ../tz2.truncoct.nc name TZ2
rms Res2_11_traj reftraj TZ2 :2-11 out rmsd.refcoords.dat
EOF
RunCpptraj "RMSD Tests."
DoTest rmsd.dat.save rmsd.dat
DoTest rmsd.mass.dat.save rmsd.mass.dat
DoTest rmsd.reftraj.dat.save rmsd.reftraj.dat
DoTest rmsd.reftraj.dat.save rmsd.refcoords.dat

# Test RMS rotate/norotate, generation of rotation matrices
TOP=""
INPUT="-i rms.in"
cat > rms.in <<EOF
parm ../tz2.parm7 [NOWAT] 
reference ../tz2.nc parm [NOWAT] 1 [first] 
parm ../tz2.truncoct.parm7 [WAT]
trajin ../tz2.truncoct.nc parm [WAT]
strip :WAT
rms NOROT ref [first] norotate @CA
outtraj tz2.norotate.crd parm [WAT]
rms ROT ref [first] out rms.dat @CA savematrices
create rmatrices.dat ROT[RM]
outtraj tz2.rotate.crd parm [WAT]
EOF
RunCpptraj "RMS coordinate rotation/rotation matrices test."
DoTest tz2.norotate.crd.save tz2.norotate.crd
DoTest tz2.rotate.crd.save tz2.rotate.crd
DoTest rmatrices.dat.save rmatrices.dat

EndTest

exit 0
