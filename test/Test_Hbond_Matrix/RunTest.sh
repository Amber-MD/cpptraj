#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nbb.dat hbavg.dat uumatrix.dat uu.dat uumatrix.gnu \
           uumatrix.nonorm.dat uumatrix.normframes.dat

INPUT="-i hbond.in"

# Solute-solute by residue matrix
TestUUresMatrix() {
  UNITNAME='Solute-solute Hbond matrix test'
  CheckFor netcdf
  if [ $? -eq 0 ] ; then
    cat > hbond.in <<EOF
noprogress
parm ../DPDP.parm7
trajin ../DPDP.nc

hbond BB out nbb.dat @N,H,C,O series uuseries uu.dat avgout hbavg.dat printatomnum \
  uuresmatrix uuresmatrixout uumatrix.nonorm.dat uuresmatrixnorm none

hbond HB uuresmatrix uuresmatrixout uumatrix.normframes.dat uuresmatrixnorm frames 

#create uumatrix.dat BB[UUresmat] nosquare2d
#create uumatrix.gnu BB[UUresmat]
run

EOF
    RunCpptraj "$UNITNAME"
    DoTest uumatrix.nonorm.dat.save uumatrix.nonorm.dat
    DoTest uumatrix.normframes.dat.save uumatrix.normframes.dat
  fi
}


TestUUresMatrix

EndTest
