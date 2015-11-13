#!/bin/bash

. ../MasterTest.sh

CleanFiles image.in reimage.mdcrd

TRAJ=ptraj.image.nc

INPUT="-i image.in"
cat > image.in <<EOF
parm ../dna30.parm7
trajin split.duplex.nc
autoimage firstatom
trajout reimage.mdcrd
EOF
RunCpptraj "AutoImage test"
DoTest reimage.mdcrd.save reimage.mdcrd
CheckTest
EndTest

exit 0 
