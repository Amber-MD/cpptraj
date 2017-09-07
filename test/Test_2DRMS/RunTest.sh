#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd1.dat rmsd2.dat ref.nc rmsd.mass.dat dme.dat trp.dat nofit.dat

RequiresNetcdf "2D RMS tests"
TOP="../tz2.parm7"
CRD="../tz2.nc"
INPUT="rms.in"

# Test 1 - 2drms
cat > rms.in <<EOF
noprogress
trajin $CRD 1 10
2drms crd1 :3-7 rmsout rmsd2.dat
EOF
RunCpptraj "2D RMSD Test."
DoTest rmsd.dat.save rmsd2.dat

# Test 2 - 2drms, mass-weighted
cat > rms.in <<EOF
noprogress
trajin $CRD 1 10 
2drms crd1 :3-7 rmsout rmsd.mass.dat mass
EOF
RunCpptraj "2D RMSD Test, mass-weighted."
DoTest rmsd.mass.dat.save rmsd.mass.dat

# Test 3 - 2drms to reference traj
TESTNAME='2D RMSD Test with reference trajectory'
CheckPnetcdf "$TESTNAME"
if [ $? -ne 0 ] ; then
  SkipCheck "$TESTNAME"
else
  cat > rms.in <<EOF
trajin $CRD 1 10
trajout ref.nc netcdf
createcrd crd1
run
2drms crdset crd1 :3-7 rmsout rmsd1.dat reftraj ref.nc
runanalysis
EOF
  RunCpptraj "$TESTNAME"
  DoTest rmsd.dat.save rmsd1.dat
fi

# Test 4 - DME
cat > rms.in <<EOF
trajin $CRD 1 10
2drms crd1 :3-7 rmsout dme.dat dme
EOF
RunCpptraj "2D DME Test."
DoTest dme.dat.save dme.dat

# Test 5 - Reference Mask
cat > rms.in <<EOF
trajin $CRD 1 10
#reference $CRD 1
#rms reference :2 :11 out trp1.dat
#trajout trp10.nc
2drms :2 :11 out trp.dat
EOF
RunCpptraj "2D RMSD with reference mask test."
DoTest trp.dat.save trp.dat

# Test 6 - No-fit RMS
cat > rms.in <<EOF
trajin $CRD 1 10
rms first :2-12@CA
2drms crd1 :2 nofit out nofit.dat
EOF
RunCpptraj "2D RMSD test, no fitting."
DoTest nofit.dat.save nofit.dat

EndTest

exit 0
