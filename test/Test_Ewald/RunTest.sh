#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ene.dat ene1.dat long.dat strip.dat directsum.0 ewald.dat \
           ew_tz2.dat ew_tz2_10.dat tz2_ortho.dat directsum.0 run9.dat

INPUT="-i ene.in"

Direct() {
  MaxThreads 1 "Direct sum test"
  if [ "$?" -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm nacl.box.parm7
trajin nacl.box.rst7
energy out directsum.0 etype directsum npoints 10 
EOF
    RunCpptraj "Direct sum test"
    DoTest directsum.0.save directsum.0
  fi
}

NaCl() {
  MaxThreads 1 "Ewald test (NaCl crystal)"
  if [ "$?" -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm nacl.box.parm7
trajin nacl.box.rst7
energy out ewald.dat etype ewald cut 5.6 dsumtol 0.0000001 \
                                 rsumtol 0.0000001 skinnb 0.01
EOF
    RunCpptraj "Ewald test (NaCl crystal)"
    DoTest ewald.dat.save ewald.dat
  fi
}

Trpzip() {
  MaxThreads 1 "Ewald test (trunc. oct)"
  if [ "$?" -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
energy out ew_tz2.dat etype ewald skinnb 0.01
EOF
    RunCpptraj "Ewald test (trunc. oct)"
    DoTest ew_tz2.dat.save ew_tz2.dat
  fi
}

Tz2_10() {
  MaxThreads 10 "Ewald test (trunc. oct), 10 frames"
  if [ "$?" -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
energy out ew_tz2_10.dat etype ewald skinnb 0.01
EOF
    RunCpptraj "Ewald test (trunc. oct), 10 frames"
    DoTest ew_tz2_10.dat.save ew_tz2_10.dat
  fi
}

Ortho() {
  MaxThreads 10 "Ewald test (Ortho), 10 frames"
  if [ "$?" -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
energy out tz2_ortho.dat etype ewald skinnb 0.01
EOF
    RunCpptraj "Ewald test (Ortho), 10 frames"
    DoTest tz2_ortho.dat.save tz2_ortho.dat
  fi
}

Direct
NaCl
Trpzip
Tz2_10
Ortho

EndTest
exit 0 
