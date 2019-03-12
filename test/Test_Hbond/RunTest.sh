#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nhb.dat avghb.dat solvhb.dat solvavg.dat solutehb2.agr \
           nbb.dat hbavg.dat solutehb.agr lifehb.gnu avg.lifehb.gnu max.lifehb.gnu \
           crv.lifehb.gnu hb?.dat hbond.mol.dat mol.avg.dat \
           ud.dat uh.dat ua.dat \
           bridgeintermol.dat avgbridgeintermol.dat noacut.dat

INPUT="-i hbond.in"

# Solute-solute, series output, lifetime analysis
TestUU() {
  UNITNAME='Solute Hbond test'
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond HB out nhb.dat avgout avghb.dat
hbond BB out nbb.dat @N,H,C,O series avgout hbavg.dat printatomnum
run
write solutehb.agr BB[solutehb] sort
runanalysis lifetime BB[solutehb] out lifehb.gnu window 10
EOF
    RunCpptraj "$UNITNAME"
    DoTest nhb.dat.save nhb.dat 
    DoTest avghb.dat.save avghb.dat
    DoTest nbb.dat.save nbb.dat
    DoTest hbavg.dat.save hbavg.dat
    DoTest solutehb.agr.save solutehb.agr
    DoTest lifehb.gnu.save lifehb.gnu
    DoTest avg.lifehb.gnu.save avg.lifehb.gnu
    DoTest max.lifehb.gnu.save max.lifehb.gnu
    DoTest crv.lifehb.gnu.save crv.lifehb.gnu
  fi
}

# Solute-solute, disk cache
TestUUcache() {
  UNITNAME='Solute Hbond test with disk caching'
  CheckFor netcdf notparallel
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc
usediskcache on
hbond BB @N,H,C,O series
run
write solutehb2.agr BB[solutehb] sort
EOF
    RunCpptraj "$UNITNAME"
    DoTest solutehb.agr.save solutehb2.agr
  fi
}

# Solute-Solvent test
TestUV() {
  UNITNAME='Solute-solvent hbond test'
  CheckFor netcdf maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
hbond hb out solvhb.dat avgout solvavg.dat :1-13 solventacceptor :WAT@O solventdonor :WAT \
      solvout solvavg.dat bridgeout solvavg.dat
EOF
    RunCpptraj "$UNITNAME"
    DoTest solvhb.dat.save solvhb.dat
    DoTest solvavg.dat.save solvavg.dat
  fi
}

# Imaged hbond test
TestImage() {
  UNITNAME='Hbond with imaging'
  CheckFor maxthreads 1
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
parm strip.4lztSc_nowat.parm7
trajin strip.4lztSc.rst7
hbond donormask :3@ND2 acceptormask :1@O avgout hb1.dat image
hbond donormask :2@ND2 acceptormask :4@O avgout hb2.dat image
EOF
    RunCpptraj "$UNITNAME"
    DoTest hb1.dat.save hb1.dat
    DoTest hb2.dat.save hb2.dat
  fi
}

# Nointramol test
TestNointramol() {
  UNITNAME='Hbond, no intramolecular hydrogen bonds test'
  CheckFor netcdf maxthreads 10
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7
trajin ../FtuFabI.NAD.TCL.nc
hbond IntraMol out hbond.mol.dat
hbond NoIntraMol out hbond.mol.dat nointramol avgout mol.avg.dat
EOF
    RunCpptraj "$UNITNAME"
    DoTest hbond.mol.dat.save hbond.mol.dat
    DoTest mol.avg.dat.save mol.avg.dat
  fi
}

# Solute specified donor mask.
SpecifiedSoluteMask() {
  UNITNAME='Hbond, specified solute masks'
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond UDmask donormask    @N=                        avgout ud.dat
hbond UHmask donormask    @N&!:PRO,NHE donorhmask @H avgout uh.dat
hbond UAmask acceptormask @O                         avgout ua.dat
EOF
    RunCpptraj "$UNITNAME"
    DoTest ud.dat.save ud.dat
    DoTest uh.dat.save uh.dat
    DoTest ua.dat.save ua.dat
  fi
}

NoAngleCut() {
  UNITNAME='Hbond, no angle cutoff test'
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond BB donormask :2-6@N acceptormask :9-13@O angle -1.0 dist 2.9 avgout noacut.dat
EOF
    RunCpptraj "$UNITNAME"
    DoTest noacut.dat.save noacut.dat
  fi
}

BridgeIntermol() {
  UNITNAME='Hbond, Bridge nointramol test'
  CheckFor netcdf maxthreads 2
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
parm ../dna30.parm7
trajin ../Test_AutoImage/split.duplex.nc
hbond hb out bridgeintermol.dat avgout avgbridgeintermol.dat :1-60 \
      solventacceptor :WAT@O solventdonor :WAT nointramol image
EOF
    RunCpptraj "$UNITNAME"
    DoTest bridgeintermol.dat.save bridgeintermol.dat
    DoTest avgbridgeintermol.dat.save avgbridgeintermol.dat
  fi
}

TestUU
TestUUcache
TestUV
TestNointramol
SpecifiedSoluteMask
NoAngleCut
BridgeIntermol
TestImage
EndTest

exit 0
