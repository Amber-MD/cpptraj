#!/bin/bash

. ../MasterTest.sh

CleanFiles change.in ala3.mod.pdb

INPUT='-i change.in'

cat > change.in <<EOF
parm ../Test_Charmm/ala3.psf
trajin ../Test_Charmm/ala3.dcd 1 1
change parmindex 0 resname from :1 to NALA
change parmindex 0 resname from :3 to CALA
change parmindex 0 atomname from @HN to H
trajout ala3.mod.pdb
EOF
RunCpptraj "Change command test"
DoTest ala3.mod.pdb.save ala3.mod.pdb

EndTest
exit 0
