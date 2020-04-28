#!/bin/bash

. ../MasterTest.sh

TESTNAME='Emin tests'

CleanFiles emin.in cpptraj.ene.dat cpptraj.emin.nc

INPUT='-i emin.in'

UNITNAME='Basic energy minimization test.'

cat > emin.in <<EOF
parm O2mol.parm7
loadcrd O2mol.rst7 name O2mol

emin crdset O2mol trajoutname cpptraj.emin.nc nsteps 100 out cpptraj.ene.dat
EOF
RunCpptraj "$UNITNAME"
DoTest cpptraj.ene.dat.save cpptraj.ene.dat

EndTest
