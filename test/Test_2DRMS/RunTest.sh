#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd1.dat rmsd2.dat ref.nc rmsd.mass.dat dme.dat trp.dat

CheckNetcdf
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
CheckTest

# Test 2 - 2drms, mass-weighted
cat > rms.in <<EOF
noprogress
trajin $CRD 1 10 
2drms crd1 :3-7 rmsout rmsd.mass.dat mass
EOF
RunCpptraj "2D RMSD Test, mass-weighted."
DoTest rmsd.mass.dat.save rmsd.mass.dat
CheckTest

# Test 3 - 2drms to reference traj
cat > rms.in <<EOF
trajin $CRD 1 10
trajout ref.nc netcdf
createcrd crd1
run
2drms crdset crd1 :3-7 rmsout rmsd1.dat reftraj ref.nc
runanalysis
EOF
RunCpptraj "2D RMSD Test with reference trajectory."
DoTest rmsd.dat.save rmsd1.dat

# Test 4 - DME
cat > rms.in <<EOF
trajin $CRD 1 10
2drms crd1 :3-7 rmsout dme.dat dme
EOF
RunCpptraj "2D DME Test."
DoTest dme.dat.save dme.dat
CheckTest

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

EndTest

exit 0
