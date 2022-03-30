#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles hbond.in nbb.dat hbavg.dat uumatrix.dat uu.dat uumatrix.gnu \
           uumatrix.nonorm.dat

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
#hbond HB out nhb.dat avgout avghb.dat
debug actions 1
hbond BB out nbb.dat @N,H,C,O series uuseries uu.dat avgout hbavg.dat printatomnum \
  uuresmatrix uuresmatrixout uumatrix.nonorm.dat uuresmatrixnorm none 
#create uumatrix.dat BB[UUresmat] nosquare2d
#create uumatrix.gnu BB[UUresmat]
run

EOF
    RunCpptraj "$UNITNAME"
    DoTest uumatrix.nonorm.dat.save uumatrix.nonorm.dat
    #DoTest nhb.dat.save nhb.dat 
    #DoTest avghb.dat.save avghb.dat
    #DoTest nbb.dat.save nbb.dat
    #DoTest hbavg.dat.save hbavg.dat
    #DoTest solutehb.agr.save solutehb.agr
    #DoTest lifehb.gnu.save lifehb.gnu
    #DoTest avg.lifehb.gnu.save avg.lifehb.gnu
    #DoTest max.lifehb.gnu.save max.lifehb.gnu
    #DoTest crv.lifehb.gnu.save crv.lifehb.gnu
  fi
}


TestUUresMatrix

EndTest
