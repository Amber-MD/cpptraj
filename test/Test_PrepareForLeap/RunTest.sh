#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in 1qos.cpptraj.pdb leap.1qos.in \
           leap.4zzw.in 4zzw.cpptraj.pdb 1qos.cpptraj.pdb.2 \
           leap.1qos.in.2 PEG.mol2 PEG.frcmod GOL.mol2 GOL.frcmod \
           PCA.mol2 PCA.frcmod MG.mol2 MG.frcmod \
           dlparams.leap.4zzw.in dlparams.4zzw.cpptraj.pdb 

INPUT='-i cpptraj.in'

# Just in case CPPTRAJHOME is not set, see if we can find the map file
RESMAPFILE=''
if [ -z "$CPPTRAJHOME" ] ; then
  if [ -e '../../dat/Carbohydrate_PDB_Glycam_Names.txt' ] ; then
    RESMAPFILE='resmapfile ../../dat/Carbohydrate_PDB_Glycam_Names.txt'
  fi
fi

cat > cpptraj.in <<EOF
parm 1qos.pdb
loadcrd 1qos.pdb name MyCrd

prepareforleap nodlparams \
  crdset MyCrd \
  name Final \
  out leap.1qos.in \
  leapunitname m \
  pdbout 1qos.cpptraj.pdb \
  nowat noh $RESMAPFILE
EOF
RunCpptraj "Prepare PDB 1qos for LEaP"
DoTest leap.1qos.in.save leap.1qos.in
DoTest 1qos.cpptraj.pdb.save 1qos.cpptraj.pdb

cat > cpptraj.in <<EOF
parm 1qos.cpptraj.pdb.save
loadcrd 1qos.cpptraj.pdb.save name MyCrd
prepareforleap nodlparams \
  hasglycam \
  crdset MyCrd \
  name Final \
  cysmask :CYX@SG \
  out leap.1qos.in.2 \
  leapunitname m \
  pdbout 1qos.cpptraj.pdb.2 \
  nowat noh $RESMAPFILE
EOF
RunCpptraj "Prepare PDB with existing glycam names (1qos)."
DoTest 1qos.cpptraj.pdb.save 1qos.cpptraj.pdb.2

cat > cpptraj.in <<EOF
parm 4zzw.pdb
loadcrd 4zzw.pdb name MyCrd

prepareforleap nodlparams \
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

# Test remote connectivity.
UNITNAME="Prepare PDB 4zzw for LEaP, download missing parameters"
PINGCMD=`which ping`
if [ ! -z "$PINGCMD" ] ; then
  $PINGCMD -c 1 raw.githubusercontent.com
  if [ $? -ne 0 ] ; then
    echo "Warning: Could not ping raw.githubusercontent.com, no internet connectivity."
    SkipCheck "$UNITNAME"
  else
    cat > cpptraj.in <<EOF
parm 4zzw.pdb
loadcrd 4zzw.pdb name MyCrd

prepareforleap dlparams bondunknown \
  crdset MyCrd \
  name Final \
  out dlparams.leap.4zzw.in \
  leapunitname m \
  pdbout dlparams.4zzw.cpptraj.pdb \
  nowat noh keepaltloc highestocc
EOF
    RunCpptraj "$UNITNAME"
    DoTest dlparams.leap.4zzw.in.save dlparams.leap.4zzw.in
  fi
else
  echo "Warning: No 'ping' found, cannot check for internet connectivity."
  SkipCheck "$UNITNAME"
fi

EndTest
