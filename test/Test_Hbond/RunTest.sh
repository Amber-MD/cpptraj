#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nhb.dat avghb.dat solvhb.dat solvavg.dat \
           nbb.dat hbavg.dat solutehb.agr lifehb.gnu avg.lifehb.gnu max.lifehb.gnu \
           crv.lifehb.gnu hb?.dat hbond.mol.dat mol.avg.dat \
           ud.dat uh.dat ua.dat

INPUT="-i hbond.in"
CheckNetcdf

# Solute-solute, series output, lifetime analysis
TestUU() {
  cat > hbond.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond HB out nhb.dat avgout avghb.dat
hbond BB out nbb.dat @N,H,C,O series avgout hbavg.dat printatomnum
run
#write solutehb.agr BB[solutehb]
write solutehb.agr BB[solutehb]:24  BB[solutehb]:104 BB[solutehb]:414 \
                   BB[solutehb]:422 BB[solutehb]:494 BB[solutehb]:652 \
                   BB[solutehb]:184 BB[solutehb]:732 BB[solutehb]:342 \
                   BB[solutehb]:352 BB[solutehb]:264 BB[solutehb]:262 \
                   BB[solutehb]:424 BB[solutehb]:218 BB[solutehb]:66 \
                   BB[solutehb]:188 BB[solutehb]:42  BB[solutehb]:86 \
                   BB[solutehb]:794
runanalysis lifetime BB[solutehb] out lifehb.gnu window 10
EOF
  RunCpptraj "Solute Hbond test."
  DoTest nhb.dat.save nhb.dat 
  DoTest avghb.dat.save avghb.dat
  DoTest nbb.dat.save nbb.dat
  DoTest hbavg.dat.save hbavg.dat
  DoTest solutehb.agr.save solutehb.agr
  DoTest lifehb.gnu.save lifehb.gnu
  DoTest avg.lifehb.gnu.save avg.lifehb.gnu
  DoTest max.lifehb.gnu.save max.lifehb.gnu
  DoTest crv.lifehb.gnu.save crv.lifehb.gnu
}

# Solute-Solvent test
TestUV() {
  cat > hbond.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
hbond hb out solvhb.dat avgout solvavg.dat :1-13 solventacceptor :WAT@O solventdonor :WAT \
      solvout solvavg.dat bridgeout solvavg.dat
EOF
  RunCpptraj "Solvent Hbond test."
  DoTest solvhb.dat.save solvhb.dat
  DoTest solvavg.dat.save solvavg.dat
}

# Imaged hbond test
TestImage() {
  MaxThreads 1 "Hbond with imaging"
  if [[ $? -eq 0 ]] ; then
    cat > hbond.in <<EOF
parm strip.4lztSc_nowat.parm7
trajin strip.4lztSc.rst7
hbond donormask :3@ND2 acceptormask :1@O avgout hb1.dat image
hbond donormask :2@ND2 acceptormask :4@O avgout hb2.dat image
EOF
    RunCpptraj "Hbond with imaging"
    DoTest hb1.dat.save hb1.dat
    DoTest hb2.dat.save hb2.dat
  fi
}

# Nointramol test
TestNointramol() {
  cat > hbond.in <<EOF
parm ../FtuFabI.NAD.TCL.parm7
trajin ../FtuFabI.NAD.TCL.nc
hbond IntraMol out hbond.mol.dat
hbond NoIntraMol out hbond.mol.dat nointramol avgout mol.avg.dat
EOF
  RunCpptraj "Hbond, no intramolecular hydrogen bonds test."
  DoTest hbond.mol.dat.save hbond.mol.dat
  DoTest mol.avg.dat.save mol.avg.dat
}

# Solute specified donor mask.
SpecifiedSoluteMask() {
  cat > hbond.in <<EOF
parm ../DPDP.parm7
trajin ../DPDP.nc
hbond UDmask donormask    @N=                        avgout ud.dat
hbond UHmask donormask    @N&!:PRO,NHE donorhmask @H avgout uh.dat
hbond UAmask acceptormask @O                         avgout ua.dat
EOF
  RunCpptraj "Hbond, specified solute masks."
  DoTest ud.dat.save ud.dat
  DoTest uh.dat.save uh.dat
  DoTest ua.dat.save ua.dat
}


TestUU
TestUV
TestImage
TestNointramol
SpecifiedSoluteMask
EndTest

exit 0
