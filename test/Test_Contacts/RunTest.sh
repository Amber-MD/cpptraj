#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in contacts.dat byres.dat byres.dat.native

NotParallel "Contacts test."
if [[ $? -ne 0 ]] ; then
  echo ""
  exit 0
fi
TOP="../tz2.truncoct.parm7"
INPUT="ptraj.in"
cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc
reference ../tz2.truncoct.nc 1

contacts first out contacts.dat :2-12@CA

contacts reference out byres.dat byresidue
EOF
RunCpptraj "Contacts test."
DoTest contacts.dat.save contacts.dat
DoTest byres.dat.save byres.dat
DoTest byres.dat.native.save byres.dat.native
CheckTest
EndTest

exit 0
