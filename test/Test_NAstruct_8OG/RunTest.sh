#!/bin/bash

. ../MasterTest.sh

CleanFiles nastruct.in *.nastruct.dat
TESTNAME='NA test (8OG)'
Requires maxthreads 2

cat > nastruct.in <<EOF
parm nowat.8OG.parm7
trajin nowat.8OG.mdcrd 1 2
nastruct naout nastruct.dat resrange 274-305 resmap 8OG:G
strip !(:274-305)
EOF
INPUT="-i nastruct.in"
RunCpptraj "$TESTNAME"
DoTest BP.nastruct.dat.save BP.nastruct.dat
DoTest BPstep.nastruct.dat.save BPstep.nastruct.dat
DoTest Helix.nastruct.dat.save Helix.nastruct.dat

EndTest
exit 0
