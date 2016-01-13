#!/bin/bash

. ../MasterTest.sh

CleanFiles charmm.in test.ala3.pdb.? test.ala3.pdb.10 first.ala3.crd \
           test.psf test.ala3.dcd second.ala3.crd
CheckNetcdf

MaxThreads 10 "Charmm DCD tests."
if [[ $? -ne 0 ]] ; then
  echo ""
  exit 0
fi

INPUT="-i charmm.in"
cat > charmm.in <<EOF
parm ala3.psf
trajin ala3.dcd 1 10
trajout test.ala3.pdb pdb multi chainid X
EOF
RunCpptraj "CHARMM PSF/DCD test"
DoTest test.ala3.pdb.save test.ala3.pdb.1
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

# Third test; psf -> psf
cat > charmm.in <<EOF
parm ala3.psf
parminfo
parmwrite out test.psf
EOF
#RunCpptraj "CHARMM PSF Write test."

EndTest

exit 0
