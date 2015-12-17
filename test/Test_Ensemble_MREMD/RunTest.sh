#!/bin/bash

. ../MasterTest.sh

CleanFiles mremd.in Strip.sorted.crd.?

INPUT="-i mremd.in"

# Test M-REMD traj sorting
cat > mremd.in <<EOF
noprogress
parm rGACC.nowat.parm7
ensemble rGACC.nowat.001
trajout Strip.sorted.crd
EOF
RunCpptraj "M-REMD sort test."
for ((i=0; i < 8; i++)) ; do
  DoTest Strip.sorted.crd.$i.save Strip.sorted.crd.$i
done

EndTest
exit 0
