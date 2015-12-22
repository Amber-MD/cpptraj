#!/bin/bash

. ../MasterTest.sh

CleanFiles mremd.in Strip.sorted.crd.? rmsd.dat rmsd.dat.? all.rmsd.dat

INPUT="-i mremd.in"

# Test M-REMD traj sorting
TrajSort() {
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
}

# Test M-REMD traj sort, rms calc
RmsTest() {
  cat > mremd.in <<EOF
noprogress
parm rGACC.nowat.parm7
ensemble rGACC.nowat.001
rms R1-4NoH first :1-4&!@H= mass out rmsd.dat
EOF
  RunCpptraj "M-REMD Rms calc test."
  if [[ -z $DO_PARALLEL ]] ; then
    DoTest rmsd.dat.save rmsd.dat
  else
    cat rmsd.dat.? > all.rmsd.dat
    DoTest all.rmsd.dat.save all.rmsd.dat
  fi
}

TrajSort
RmsTest

EndTest
exit 0
