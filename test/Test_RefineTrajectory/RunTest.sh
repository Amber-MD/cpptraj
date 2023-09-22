#!/bin/bash

. ../MasterTest.sh

CleanFiles cpptraj.in rms.dat refinedcoords.crd normcoords.crd normcoords2.crd \
           trimcoords.crd trimcoords2.crd

INPUT='-i cpptraj.in'

TESTNAME='Trajectory refinement tests'

Requires netcdf

# --------------------------------------
RefineWithLoop() {
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
reference ../tz2.nc 1 [FIRST]

# Initial average
rms first !@H=
average crdset R0
run

crdaction R0 rms Rinit ref [FIRST] out rms.dat !@H=
# Loop over references
for i=0;i<9;i++
  j = \$i + 1
  rms ref R\$i !@H=
  average crdset R\$j
  trajout looprefined.crd
  run
  crdaction R\$j rms R\$i.\$j ref R\$i out rms.dat !@H=
done
#crdout R\$j looprefined.crd
list
EOF
  RunCpptraj "$TESTNAME, RMS refinement using loop"
  DoTest refinedcoords.crd.save looprefined.crd -a 0.002
}

# --------------------------------------
RmsRefine() {
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd

crdtransform MyCrd rmsrefine mask !@H=
crdout MyCrd refinedcoords.crd 
EOF
  RunCpptraj "$TESTNAME, using crdtransform rmsrefine"
  DoTest refinedcoords.crd.save refinedcoords.crd
}

# --------------------------------------
NormCoords() {
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd

crdtransform MyCrd normcoords
crdout MyCrd normcoords.crd
EOF
  RunCpptraj "$TESTNAME, using crdtransform normcoords"
  DoTest normcoords.crd.save normcoords.crd

  cat > cpptraj.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd

crdtransform MyCrd normcoords name NewCrd
crdout NewCrd normcoords2.crd
EOF
  RunCpptraj "$TESTNAME, using crdtransform normcoords to new set"
  DoTest normcoords.crd.save normcoords2.crd
}

# --------------------------------------
Trim() {
  cat > cpptraj.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd

crdtransform MyCrd trim cutoff 0.1
crdout MyCrd trimcoords.crd
list
EOF
  RunCpptraj "$TESTNAME, using crdtransform trim"
  DoTest trimcoords.crd.save trimcoords.crd

  cat > cpptraj.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.nc name MyCrd

crdtransform MyCrd trim cutoff 0.1 criterion medoid
crdout MyCrd trimcoords2.crd
list
EOF
  RunCpptraj "$TESTNAME, using crdtransform trim criterion medoid"
  DoTest trimcoords2.crd.save trimcoords2.crd
}

# --------------------------------------
RefineWithLoop
RmsRefine
NormCoords
Trim

EndTest

