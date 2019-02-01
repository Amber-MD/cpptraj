#!/bin/bash

. ../MasterTest.sh

CleanFiles matrix.dat hd.in hd.dat

INPUT='-i hd.in'

cat > matrix.dat <<EOF
1.0000       2.2361       3.0000       4.1231
2.2361       1.0000       2.2361       5.0000
3.0000       2.2361       1.0000       4.1231
2.2361       3.0000       2.2361       3.0000
EOF

cat > hd.in <<EOF
readdata matrix.dat read2d name Matrix
runanalysis hausdorff Matrix out hd.dat name HD outab hd.dat outba hd.dat
EOF
RunCpptraj "Simple Hausdorff distance test."
DoTest hd.dat.save hd.dat

EndTest
exit 0
