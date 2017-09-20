#!/bin/bash

. ../MasterTest.sh

CleanFiles mremd.in Strip.sorted.crd.? rmsd.dat rmsd.dat.? all.rmsd.dat \
           nhbond.dat nhbond.dat.? all.nhbond.dat \
           hbavg.dat hbavg.dat.? all.hbavg.dat \
           Outtraj.crd Outtraj.crd.0 Outtraj.crd.1 avg.rst7.? \
           RA.dat RA.dat.? all.RA.dat melt.dat RmsToRep1.dat \
           RmsToRep1.dat.?

INPUT="-i mremd.in"
TESTNAME='M-REMD ensemble tests'
Requires netcdf nthreads 8

# Test M-REMD traj sorting
TrajSort() {
  cat > mremd.in <<EOF
noprogress
parm rGACC.nowat.parm7
ensemblesize 8
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
ensemblesize 8
ensemble rGACC.nowat.001
hbond HB :1-4 solventdonor :Na+ solventacceptor :Na+ \
      out nhbond.dat avgout hbavg.dat
rms R1-4NoH first :1-4&!@H= mass out rmsd.dat
average avg.rst7 :1-4
run
EOF
  RunCpptraj "M-REMD actions test."
  if [ -z "$DO_PARALLEL" ] ; then
    DoTest rmsd.dat.save rmsd.dat
    DoTest nhbond.dat.save nhbond.dat
    DoTest hbavg.dat.save hbavg.dat
  else
    cat rmsd.dat.? > all.rmsd.dat     && DoTest all.rmsd.dat.save all.rmsd.dat
    cat nhbond.dat.? > all.nhbond.dat && DoTest all.nhbond.dat.save all.nhbond.dat
    cat hbavg.dat.? > all.hbavg.dat   && DoTest all.hbavg.dat.save all.hbavg.dat
  fi
  DoTest avg.rst7.2.save avg.rst7.2
  DoTest avg.rst7.5.save avg.rst7.5
}

# Test M-REMD traj sort, actions + analysis
AnalysisTest() {
  cat > mremd.in <<EOF
noprogress
parm rGACC.nowat.parm7
reference rGACC.nowat.001
ensemblesize 8
ensemble rGACC.nowat.001
rms Ref reference out RmsToRep1.dat
EOF
  RunCpptraj "M-REMD, generate RMS data"
  if [ -z "$DO_PARALLEL" ] ; then
    READDATA='readdata RmsToRep1.dat name Ref'
  else
    READDATA='readdata RmsToRep1.dat.? name Ref separate'
  fi
  cat > mremd.in <<EOF
$READDATA
runanalysis meltcurve Ref* out melt.dat cut 5.0 name ToRep1
EOF
  RunCpptraj "M-REMD, meltcurve analysis"
  DoTest melt.dat.save melt.dat
}

# Test M-REMD process, no sort, running average (tests preload in parallel)
RunAvgTest() {
  cat > mremd.in <<EOF
parm rGACC.nowat.parm7
ensemblesize 8
ensemble rGACC.nowat.001 nosort
runavg window 3
rms RA first :1-4&!@H= out RA.dat
EOF
  RunCpptraj "M-REMD no sort, running average test"
  if [ -z "$DO_PARALLEL" ] ; then
    DoTest RA.dat.save RA.dat
  else
    cat RA.dat.? > all.RA.dat && DoTest all.RA.dat.save all.RA.dat
  fi
}

# Test ensemble mode with outtraj output
OuttrajTest() {
  cat > mremd.in <<EOF
noprogress
parm rGACC.nowat.parm7
ensemblesize 8
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

TrajSort
ActionsTest
OuttrajTest
RunAvgTest
AnalysisTest

EndTest
exit 0
