#!/bin/bash

. ../MasterTest.sh

CleanFiles pw.in pair.dat energy.dat avg.dat map.elec.gnu map.vdw.gnu

TESTNAME='Pairwise test'
Requires netcdf

INPUT='-i pw.in'

cat > pw.in <<EOF
parm ../Test_SymmRmsd/TYR.parm7
trajin ../Test_SymmRmsd/TYR.nc
pairwise TYR out pair.dat eout energy.dat avgout avg.dat \
         vmapout map.vdw.gnu emapout map.elec.gnu \
         printmode or cuteelec 2.0 cutevdw 2.0
EOF
RunCpptraj "$TESTNAME"
DoTest pair.dat.save pair.dat
DoTest energy.dat.save energy.dat
DoTest avg.dat.save avg.dat
DoTest map.vdw.gnu.save map.vdw.gnu
DoTest map.elec.gnu.save map.elec.gnu

EndTest
exit 0
