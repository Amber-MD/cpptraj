#!/bin/bash

. ../MasterTest.sh

INPUT="diffusion.in"
TOP=../tz2.ortho.parm7

CleanFiles $INPUT diff_?.xmgr diff.dat diff.1.dat diff.2.dat diff.3.dat nw.dat \
           WAT_O.agr DC.dat
CheckNetcdf
# Basic ptraj diffusion test
# creates <prefix>_X.xmgr, X = {a,r,x,y,z}
Test_diffusion_oldSyntax() {
  cat > $INPUT <<EOF
trajin ../tz2.ortho.nc
diffusion :1-13 1.0 diff
EOF
  RunCpptraj "Diffusion test, old syntax."
  DoTest diff_a.xmgr.save diff_a.xmgr
  DoTest diff_r.xmgr.save diff_r.xmgr
  DoTest diff_x.xmgr.save diff_x.xmgr
  DoTest diff_y.xmgr.save diff_y.xmgr
  DoTest diff_z.xmgr.save diff_z.xmgr
}

Test_diffusion_newSyntax() {
cat > $INPUT <<EOF
trajin ../tz2.ortho.nc
diffusion :WAT@O out WAT_O.agr WAT_O diffout DC.dat
EOF
  RunCpptraj "Diffusion test, new syntax."
  DoTest WAT_O.agr.save WAT_O.agr -r 0.0007
  DoTest DC.dat.save DC.dat -r 0.004
}

# ------------------------------------------------------------------------------
# STFC diffusion tests
Test_stfc_diffusion() {
  CMD="stfcdiffusion"
  # Basic test
  cat > $INPUT <<EOF
trajin ../tz2.ortho.nc
$CMD mask :1-13 out diff.dat
EOF
  RunCpptraj "STFC Diffusion Test"
  DoTest diff.dat.save diff.dat
  # Test with individual distances
  cat > $INPUT <<EOF
trajin ../tz2.ortho.nc
$CMD mask :1-13 out diff.1.dat distances
EOF
  RunCpptraj "STFC Diffusion Test with individual distances"
  DoTest diff.1.dat.save diff.1.dat
  # Test with COM
  cat > $INPUT <<EOF
trajin ../tz2.ortho.nc
$CMD mask :1-13 out diff.2.dat com
EOF
  RunCpptraj "STFC Diffusion Test with COM"
  DoTest diff.2.dat.save diff.2.dat
  # Test with second mask
  cat > $INPUT <<EOF
trajin ../tz2.ortho.nc 
$CMD mask :1-13 out diff.3.dat mask2 :WAT
EOF
  RunCpptraj "STFC Diffusion Test with second mask"
  DoTest diff.3.dat.save diff.3.dat
  DoTest nw.dat.save nw.dat
}

Test_diffusion_oldSyntax
Test_diffusion_newSyntax
NotParallel "STFC diffusion tests."
if [[ $? -eq 0 ]] ; then
  Test_stfc_diffusion
fi

EndTest

exit 0
