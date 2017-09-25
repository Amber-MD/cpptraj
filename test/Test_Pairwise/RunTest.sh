#!/bin/bash

. ../MasterTest.sh

CleanFiles pw.in pair.dat energy.dat avg.dat map.elec.gnu map.vdw.gnu \
           ref.ene.pdb ref.cut.mol2.*.mol2.?

TESTNAME='Pairwise test'
Requires netcdf maxthreads 1

INPUT='-i pw.in'

cat > pw.in <<EOF
parm ../Test_SymmRmsd/TYR.parm7
trajin ../Test_SymmRmsd/TYR.nc
pairwise TYR out pair.dat eout energy.dat avgout avg.dat \
         vmapout map.vdw.gnu emapout map.elec.gnu \
         printmode or cuteelec 2.0 cutevdw 2.0
run

clear trajin
reference ../Test_SymmRmsd/TYR.nc 1
trajin ../Test_SymmRmsd/TYR.nc 2 2
pairwise REF reference cutout ref.cut.mol2 pdbout ref.ene.pdb scalepdbe
EOF
RunCpptraj "$TESTNAME"
DoTest pair.dat.save pair.dat
DoTest energy.dat.save energy.dat
DoTest avg.dat.save avg.dat
DoTest map.vdw.gnu.save map.vdw.gnu
DoTest map.elec.gnu.save map.elec.gnu
DoTest ref.cut.mol2.eelec.mol2.1.save ref.cut.mol2.eelec.mol2.1
DoTest ref.cut.mol2.evdw.mol2.1.save ref.cut.mol2.evdw.mol2.1
DoTest ref.ene.pdb.save ref.ene.pdb

EndTest
exit 0
