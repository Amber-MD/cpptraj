#!/bin/bash

. ../MasterTest.sh

CleanFiles matrix.in mtest.dat.save mtest.*.dat dihcovar.dat modes.dihcovar.dat \
           dih.project.dat dih.project.agr dihedrals.dat

TESTNAME='Dihedral Covariance Matrix Test'
RequiresMathlib "$TESTNAME" 

INPUT="-i matrix.in"
if [ -z "$DO_PARALLEL" ] ; then
  # Serial version. Can do everything in one run.
  cat > matrix.in <<EOF
parm ../Test_Matrix/1rrb_vac.prmtop
trajin ../Test_Matrix/1rrb_vac.mdcrd

multidihedral BB phi psi resrange 2
run
matrix dihcovar dihedrals BB[*] out dihcovar.dat name DIH
precision dihcovar.dat 12 6 
diagmatrix DIH vecs 4 out modes.dihcovar.dat name DIHMODES
run
projection evecs DIHMODES out dih.project.dat beg 1 end 4 dihedrals BB[*]
#create dih.project.agr Proj_00004 
EOF
  RunCpptraj "Dihedral Covariance Matrix Test"
else
  # Parallel version. Need two runs.
  cat > matrix.in <<EOF
parm ../Test_Matrix/1rrb_vac.prmtop
trajin ../Test_Matrix/1rrb_vac.mdcrd
multidihedral BB phi psi resrange 2 out dihedrals.dat
run
matrix dihcovar dihedrals BB[*] out dihcovar.dat name DIH
precision dihcovar.dat 12 6 
diagmatrix DIH vecs 4 out modes.dihcovar.dat name DIHMODES
run
EOF
  RunCpptraj "Dihedral Covariance Matrix Test (parallel, create matrix)"
  cat > matrix.in <<EOF
parm ../Test_Matrix/1rrb_vac.prmtop
trajin ../Test_Matrix/1rrb_vac.mdcrd
readdata dihedrals.dat name BB 
dataset mode torsion BB
projection evecs modes.dihcovar.dat out dih.project.dat beg 1 end 4 dihedrals BB
EOF
  RunCpptraj "Dihedral Covariance Matrix Test (parallel, projection)"
fi
DoTest dihcovar.dat.save dihcovar.dat
DoTest modes.dihcovar.dat.save modes.dihcovar.dat
DoTest dih.project.dat.save dih.project.dat

EndTest
  
exit 0
