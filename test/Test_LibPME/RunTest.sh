#!/bin/bash

. ../MasterTest.sh

CleanFiles ene.in ewald.dat debug.nacl.dat nacl.dat debug.tz2n.dat tz2n.dat \
           debug.tz2o.dat tz2o.dat debug.mtz2o.dat mtz2o.dat pme.nacl.dat \
           long_tz2n.dat
INPUT="-i ene.in"
TESTNAME='Particle mesh Ewald tests'
Requires libpme maxthreads 10
# Set to 1 for debugging purposes
PMEDEBUG=0
if [ $PMEDEBUG -eq 0 ] ; then
  ECMD='#energy'
  PREFIX=''
else
  ECMD='energy'
  PREFIX='debug.'
fi

Simple() {
  UNITNAME='Particle mesh Ewald test (simple)'
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm test.mol2
trajin test.mol2
box x 20 y 20 z 20 alpha 90 beta 90 gamma 90
energy out ewald.dat etype pme cut 5.6 dsumtol 0.0000001 skinnb 0.01
#vector UX ucellx
#vector UY ucelly
#vector UZ ucellz
#run
#writedata ucell.mol2 vectraj trajfmt mol2 UX UY UZ
EOF
    RunCpptraj "$UNITNAME"
  fi
}

NaCl() {
  UNITNAME='Particle mesh Ewald test (NaCl crystal)'
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    TFILE="$PREFIX"nacl.dat
    cat > ene.in <<EOF
noprogress
parm ../Test_Ewald/nacl.box.parm7
trajin ../Test_Ewald/nacl.box.rst7
debug actions $PMEDEBUG
$ECMD Reg nonbond out $TFILE etype ewald cut 5.6 dsumtol 0.0000001 rsumtol 0.000000001 skinnb 0.01 mlimits 12,12,12

energy Pme nonbond out $TFILE etype pme cut 5.6 dsumtol 0.0000001 skinnb 0.01 nfft 32,32,32
EOF
    RunCpptraj "$UNITNAME"
    if [ $PMEDEBUG -gt 0 ] ; then
      grep "DEBUG: Eself" test.out > pme.nacl.dat
      DoTest pme.nacl.dat.save pme.nacl.dat
    fi
    DoTest "$TFILE".save "$TFILE"
  fi
}

TrpzipNonortho() {
  UNITNAME='Particle mesh Ewald test (trunc. oct)'
  CheckFor netcdf maxthreads 1
  if [ $? -eq 0 ] ; then
    TFILE="$PREFIX"tz2n.dat
    cat > ene.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1
#debug actions 1
$ECMD Reg nonbond out $TFILE etype ewald skinnb 0.01 \
       cut 8.0 dsumtol 0.0000001 rsumtol 0.000000001
energy Pme nonbond out $TFILE etype pme   skinnb 0.01 order 6 \
       cut 8.0 dsumtol 0.0000001 nfft 96,90,90
precision $TFILE 20 10
EOF
    RunCpptraj "$UNITNAME"
    DoTest tz2n.dat.save tz2n.dat -a 0.000001
  fi
}

TrpzipOrtho() {
  UNITNAME='Particle mesh Ewald test (ortho)'
  CheckFor netcdf maxthreads 1
  if [ $? -eq 0 ] ; then
    TFILE="$PREFIX"tz2o.dat
    cat > ene.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
#debug actions 1
$ECMD Reg nonbond out $TFILE etype ewald skinnb 0.01 \
       cut 8.0 dsumtol 0.0000001 rsumtol 0.000000001
energy Pme nonbond out $TFILE etype pme   skinnb 0.01 order 6 \
       cut 8.0 dsumtol 0.0000001 nfft 72,90,72
precision $TFILE 20 10
EOF
    RunCpptraj "$UNITNAME"
    DoTest "$TFILE".save "$TFILE" -a 0.000001
  fi
}

MaskTz2Ortho() {
  UNITNAME='Particle mesh Ewald test (ortho, with mask)'
  CheckFor netcdf maxthreads 1
  if [ $? -eq 0 ] ; then
    TFILE="$PREFIX"mtz2o.dat
    cat > ene.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 1
#debug actions 1
$ECMD Reg nonbond out $TFILE etype ewald skinnb 0.01 !:WAT \
       cut 8.0 dsumtol 0.0000001 rsumtol 0.000000001
energy Pme nonbond out $TFILE etype pme   skinnb 0.01 order 6 !:WAT \
       cut 8.0 dsumtol 0.0000001 nfft 72,90,72
precision $TFILE 20 10
EOF
    RunCpptraj "$UNITNAME"
    DoTest "$TFILE".save "$TFILE" -a 0.0000001
  fi
}

Tz2_Nonortho_10() {
  UNITNAME='PME test (trunc. oct), 10 frames'
  CheckFor netcdf long
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
energy Pme nonbond out long_tz2n.dat etype pme skinnb 2.0 cut 8.0 \
  dsumtol 0.0000001 nfft 72,90,72
EOF
    RunCpptraj "$UNITNAME"
    DoTest long_tz2n.dat.save long_tz2n.dat
  fi
}

Tz2_Ortho_10() {
  UNITNAME='Ewald test (ortho), 10 frames'
  CheckFor netcdf long
  if [ $? -eq 0 ] ; then
    cat > ene.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
energy out tz2_ortho.dat etype ewald skinnb 0.01
EOF
    RunCpptraj "Ewald test (ortho), 10 frames"
    DoTest tz2_ortho.dat.save tz2_ortho.dat
  fi
}

#Simple
NaCl
TrpzipNonortho
TrpzipOrtho
MaskTz2Ortho
Tz2_Nonortho_10
#Tz2_Ortho_10

EndTest
exit 0 
