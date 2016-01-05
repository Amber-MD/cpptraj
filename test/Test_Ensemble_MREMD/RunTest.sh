#!/bin/bash

. ../MasterTest.sh

CleanFiles mremd.in Strip.sorted.crd.? rmsd.dat rmsd.dat.? all.rmsd.dat \
           nhbond.dat nhbond.dat.? all.nhbond.dat \
           hbavg.dat hbavg.dat.? all.hbavg.dat \
           Outtraj.crd Outtraj.crd.0 Outtraj.crd.1

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

# Test M-REMD traj sort, actions 
ActionsTest() {
  cat > mremd.in <<EOF
noprogress
parm rGACC.nowat.parm7
ensemble rGACC.nowat.001
hbond HB :1-4 solventdonor :Na+ solventacceptor :Na+ \
      out nhbond.dat avgout hbavg.dat
rms R1-4NoH first :1-4&!@H= mass out rmsd.dat
EOF
  RunCpptraj "M-REMD actions test."
  if [[ -z $DO_PARALLEL ]] ; then
    DoTest rmsd.dat.save rmsd.dat
    DoTest nhbond.dat.save nhbond.dat
    DoTest hbavg.dat.save hbavg.dat
  else
    cat rmsd.dat.? > all.rmsd.dat     && DoTest all.rmsd.dat.save all.rmsd.dat
    cat nhbond.dat.? > all.nhbond.dat && DoTest all.nhbond.dat.save all.nhbond.dat
    cat hbavg.dat.? > all.hbavg.dat   && DoTest all.hbavg.dat.save all.hbavg.dat
  fi
  
}

# Test ensemble mode with outtraj output
OuttrajTest() {
  cat > mremd.in <<EOF
noprogress
parm rGACC.nowat.parm7
ensemble rGACC.nowat.001 nosort
outtraj Outtraj.crd onlymembers 0,1
EOF
  RunCpptraj "M-REMD ensemble outtraj test."
  if [[ -e Outtraj.crd ]] ; then # FIXME should be unnecessary when merged back to master
    DoTest Outtraj.crd.0.save Outtraj.crd
  else
    DoTest Outtraj.crd.0.save Outtraj.crd.0
  fi
  DoTest Outtraj.crd.1.save Outtraj.crd.1
}

RequiresThreads 8
if [[ $? -eq 0 ]] ; then
  TrajSort
  ActionsTest
  OuttrajTest
else
  echo ""
  exit 0
fi

EndTest
exit 0
