#!/bin/bash

. ../MasterTest.sh

CleanFiles fix.in order.in reorder.outoforder.parm7 reorder.mdcrd

INPUT='-i fix.in'

cat > fix.in <<EOF
readinput order.in
EOF

cat > order.in <<EOF
parm outoforder.parm7
trajin min1.crd 1 10
fixatomorder outprefix reorder
trajout reorder.mdcrd
EOF
RunCpptraj "Out of order molecules test"
DoTest reorder.outoforder.parm7.save reorder.outoforder.parm7 -I %VERSION
DoTest reorder.mdcrd.save reorder.mdcrd

EndTest
exit 0
