#!/bin/bash

. ../MasterTest.sh

INPUT="diffusion.in"

CleanFiles $INPUT diff_?.xmgr diff.dat diff.1.dat diff.2.dat diff.3.dat nw.dat \
           WAT_O.agr DC.dat Nonortho.agr Nonortho.dat noimage.agr noimage.dat
CheckNetcdf

# Basic ptraj diffusion test
# creates <prefix>_X.xmgr, X = {a,r,x,y,z}
Test_diffusion_oldSyntax() {
  TOP=../tz2.ortho.parm7
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
  TOP=../tz2.ortho.parm7
  cat > $INPUT <<EOF
trajin ../tz2.ortho.nc
diffusion :WAT@O out WAT_O.agr WAT_O diffout DC.dat
EOF
  RunCpptraj "Diffusion test, new syntax."
  DoTest WAT_O.agr.save WAT_O.agr
  DoTest DC.dat.save DC.dat
}

Test_diffusion_nonOrtho() {
  TOP=../tz2.truncoct.parm7
  cat > $INPUT <<EOF
trajin ../tz2.truncoct.nc
diffusion :WAT@O out Nonortho.agr WAT_O diffout Nonortho.dat
EOF
  RunCpptraj "Diffusion test, nonorthogonal box"
  DoTest Nonortho.dat.save Nonortho.dat 
  DoTest Nonortho.agr.save Nonortho.agr
}

Test_diffusion_noImage() {
  TOP=../tz2.truncoct.parm7
  cat > $INPUT <<EOF
trajin ../tz2.truncoct.nc
diffusion :WAT@O out noimage.agr WAT_O diffout noimage.dat noimage
EOF
  RunCpptraj "Diffusion test, no imaging."
  DoTest noimage.dat.save noimage.dat
  DoTest noimage.agr.save noimage.agr
}

# ------------------------------------------------------------------------------
# STFC diffusion tests
Test_stfc_diffusion() {
  TOP=../tz2.ortho.parm7
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

Test_diffusion_noImage
MaxThreads 1 "Imaged diffusion tests"
if [[ $? -eq 0 ]] ; then
  Test_diffusion_oldSyntax
  Test_diffusion_newSyntax
  Test_diffusion_nonOrtho
fi
NotParallel "STFC diffusion tests."
if [[ $? -eq 0 ]] ; then
  Test_stfc_diffusion
fi

EndTest

exit 0
