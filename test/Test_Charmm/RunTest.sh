#!/bin/bash

. ../MasterTest.sh

CleanFiles charmm.in test.ala3.pdb.? test.ala3.pdb.10 first.ala3.crd \
           test.psf test.ala3.dcd second.ala3.crd strip.chamber.parm7 \
           run0.res_0.mol2 cpptraj.psf

TESTNAME='Charmm DCD tests'
Requires maxthreads 10

INPUT="-i charmm.in"
cat > charmm.in <<EOF
parm ala3.psf
trajin ala3.dcd 1 10
trajout test.ala3.pdb pdb multi chainid X
EOF
RunCpptraj "CHARMM PSF/DCD test"
DoTest test.ala3.pdb.save test.ala3.pdb.1

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

# CHAMBER strip test
cat > charmm.in <<EOF
parm chamber.ala3.parm7
parmstrip :2-3
parmwrite out strip.chamber.parm7
quit
EOF
RunCpptraj "CHAMBER topology read/strip test."
DoTest strip.chamber.parm7.save strip.chamber.parm7 -I %VERSION

UNITNAME='Read CHARMM restart'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > charmm.in <<EOF
parm ala3.psf
trajin run0.res_0
trajout run0.res_0.mol2
EOF
  RunCpptraj "$UNITNAME"
  DoTest run0.res_0.mol2.save run0.res_0.mol2
fi

UNITNAME='CHARMM PSF write'
cat > charmm.in <<EOF
parm ala3.psf
parmwrite out cpptraj.psf
EOF
RunCpptraj "$UNITNAME"

EndTest

exit 0
