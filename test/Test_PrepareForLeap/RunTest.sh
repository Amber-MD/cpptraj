#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in 1qos.cpptraj.pdb leap.1qos.in \
           leap.4zzw.in 4zzw.cpptraj.pdb

INPUT='-i cpptraj.in'

cat > cpptraj.in <<EOF
parm 1qos.pdb
loadcrd 1qos.pdb name MyCrd

prepareforleap \
  crdset MyCrd \
  name Final \
  out leap.1qos.in \
  leapunitname m \
  pdbout 1qos.cpptraj.pdb \
  nowat noh
EOF
RunCpptraj "Prepare PDB 1qos for LEaP"
DoTest leap.1qos.in.save leap.1qos.in
DoTest 1qos.cpptraj.pdb.save 1qos.cpptraj.pdb

cat > cpptraj.in <<EOF
parm 4zzw.pdb
loadcrd 4zzw.pdb name MyCrd

prepareforleap \
  crdset MyCrd \
  name Final \
  out leap.4zzw.in \
  leapunitname m \
  pdbout 4zzw.cpptraj.pdb \
  nowat noh keepaltloc highestocc
EOF
RunCpptraj "Prepare PDB 4zzw for LEaP"
DoTest leap.4zzw.in.save leap.4zzw.in
DoTest 4zzw.cpptraj.pdb.save 4zzw.cpptraj.pdb

EndTest
