#!/bin/bash

. ../MasterTest.sh

CleanFiles charmm.in test.ala3.pdb test.ala3.nc first.ala3.crd test.ala3.dcd second.ala3.crd
CheckNetcdf
cat > charmm.in <<EOF
parm ala3.psf
trajin ala3.dcd
trajout test.ala3.pdb pdb onlyframes 1 chainid X
trajout test.ala3.nc netcdf
EOF
INPUT="-i charmm.in"
RunCpptraj "CHARMM PSF/DCD test"
DoTest test.ala3.pdb.save test.ala3.pdb
CheckTest

# Second test: Read in 10 frames of a dcd traj, write
# both an Amber coord and dcd traj. Then read the written
# dcd traj and convert to second amber traj. Compare the
# written amber trajs
cat > charmm.in <<EOF
parm ala3.psf
trajin ala3.dcd 1 10
trajout first.ala3.crd
trajout test.ala3.dcd dcd
EOF
RunCpptraj "CHARMM DCD Write, step 1."

cat > charmm.in <<EOF
parm ala3.psf
trajin test.ala3.dcd
trajout second.ala3.crd
EOF
RunCpptraj "CHARMM DCD Write, step 2."
DoTest first.ala3.crd second.ala3.crd

EndTest

exit 0
