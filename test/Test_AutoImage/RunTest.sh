#!/bin/bash

. ../MasterTest.sh

CleanFiles image.in reimage.mdcrd
CheckNetcdf
TRAJ=ptraj.image.nc
INPUT="-i image.in"
MaxThreads 2 "AutoImage test"
if [[ $? -eq 0 ]] ; then
  cat > image.in <<EOF
parm ../dna30.parm7
trajin split.duplex.nc
autoimage firstatom
trajout reimage.mdcrd
EOF
  RunCpptraj "AutoImage test"
  DoTest reimage.mdcrd.save reimage.mdcrd
fi
EndTest

exit 0 
