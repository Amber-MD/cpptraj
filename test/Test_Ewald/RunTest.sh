#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ene.dat ene1.dat long.dat strip.dat directsum.0 ewald.dat \
           ew_tz2.dat ew_tz2_10.dat tz2_ortho.dat directsum.0 run9.dat

INPUT="-i ene.in"

Direct() {
  cat > ene.in <<EOF
noprogress
parm nacl.box.parm7
trajin nacl.box.rst7
energy out directsum.0 etype directsum npoints 10 
EOF
  RunCpptraj "Direct sum test"
  DoTest directsum.0.save directsum.0
}

NaCl() {
  cat > ene.in <<EOF
noprogress

parm nacl.box.parm7
trajin nacl.box.rst7
energy out ewald.dat etype ewald cut 5.6 dsumtol 0.0000001 \
                                 rsumtol 0.0000001 skinnb 0.01
EOF
  RunCpptraj "Ewald test (NaCl crystal)"
  DoTest ewald.dat.save ewald.dat
}

Trpzip() {
  cat > ene.in <<EOF
noprogress

parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
energy out ew_tz2.dat etype ewald skinnb 0.01
EOF
  RunCpptraj "Ewald test (trunc. oct)"
  DoTest ew_tz2.dat.save ew_tz2.dat
}

Tz2_10() {
  cat > ene.in <<EOF
noprogress

parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
energy out ew_tz2_10.dat etype ewald skinnb 0.01
EOF
  RunCpptraj "Ewald test (trunc. oct), 10 frames"
  DoTest ew_tz2_10.dat.save ew_tz2_10.dat
}

Ortho() {
  cat > ene.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
energy out tz2_ortho.dat etype ewald skinnb 0.01
EOF
  RunCpptraj "Ewald test (Ortho)"
  DoTest tz2_ortho.dat.save tz2_ortho.dat
}

FtuFabI() {
  cat > ene.in <<EOF
noprogress
parm ../ChainA_1-268_NAD_TCL-gaff.tip3p.parm7
trajin ../run9.nc 1 1
energy out run9.dat etype ewald
EOF
  RunCpptraj "Ewald test (FtuFabI)"
  DoTest run9.dat.save run9.dat
}

Direct
NaCl
Trpzip
Tz2_10
Ortho
if [ ! -z "$OPENMP" ] ; then
  FtuFabI
else
  echo ""
  echo "FtuFabI test only run for OpenMP."
  echo ""
fi

EndTest
exit 0 
