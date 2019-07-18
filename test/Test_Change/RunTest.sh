#!/bin/bash

. ../MasterTest.sh

CleanFiles change.in ala3.mod.pdb ala3.chain.pdb crdala3.chain.pdb

TESTNAME='Change command test'
Requires maxthreads 1

INPUT='-i change.in'

cat > change.in <<EOF
parm ../Test_Charmm/ala3.psf
trajin ../Test_Charmm/ala3.dcd 1 1
change parmindex 0 resname from :1 to NALA
change parmindex 0 resname from :3 to CALA
change parmindex 0 atomname from @HN to H
trajout ala3.mod.pdb
EOF
RunCpptraj "Change atom and residue name tests"
DoTest ala3.mod.pdb.save ala3.mod.pdb

cat > change.in <<EOF
parm ../Test_Charmm/ala3.psf
trajin ../Test_Charmm/ala3.dcd 1 1
change parmindex 0 chainid of * to A
trajout ala3.chain.pdb
EOF
RunCpptraj "Change chain ID test"
DoTest ala3.chain.pdb.save ala3.chain.pdb

cat > change.in <<EOF
parm ../Test_Charmm/ala3.psf
loadcrd ../Test_Charmm/ala3.dcd 1 1 name MyCrd
change crdset MyCrd chainid of * to A
crdout MyCrd crdala3.chain.pdb
EOF
RunCpptraj "Change chain ID of COORDS set test"
DoTest ala3.chain.pdb.save crdala3.chain.pdb

EndTest
exit 0
