#!/bin/bash

. ../MasterTest.sh

TESTNAME='XYZ format tests'

CleanFiles xyz.in tz2.xyz test1.crd.save test?.crd tz2.st.xyz tz2.nt.at.xyz \
           tz2.nt.xyz tz2.mt.at.xyz tz2.mt.xyz

INPUT='-i xyz.in'

UNITNAME='XYZ format write'
cat > xyz.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd 1 10
trajout test1.crd.save crd
# Single title, atom-xyz
trajout tz2.xyz
# Single title, xyz
trajout tz2.st.xyz    titletype single   ftype xyz
# No title, atom-xyz
trajout tz2.nt.at.xyz titletype none     ftype atomxyz
# No title, xyz
trajout tz2.nt.xyz    titletype none     ftype xyz
# Multi title, atom-xyz
trajout tz2.mt.at.xyz titletype perframe ftype atomxyz
# Multi title, xyz
trajout tz2.mt.xyz    titletype perframe ftype xyz
EOF
RunCpptraj "$UNITNAME"

UNITNAME='Atom XYZ format read'
N=1
for FILE in tz2.xyz tz2.st.xyz tz2.nt.at.xyz tz2.nt.xyz tz2.mt.at.xyz tz2.mt.xyz ; do
  cat > xyz.in <<EOF
parm ../tz2.parm7
trajin $FILE 
trajout test$N.crd
EOF
  RunCpptraj "$UNITNAME, test $N"
  DoTest test1.crd.save test$N.crd
  ((N++))
done

EndTest
exit 0
