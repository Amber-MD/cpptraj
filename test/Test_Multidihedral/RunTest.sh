#!/bin/bash

. ../MasterTest.sh

CleanFiles dih.in multidih.dat dihedral.dat custom.dat dihedral2.dat all.dat \
           chin.dat arg.dat cyclic.dat

TESTNAME='Multidihedral tests'
Requires maxthreads 10

INPUT="-i dih.in"
UNITNAME='Multidihedral test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc 1 10
multidihedral out multidih.dat phi psi
EOF
  for ((res = 1; res <= 12; res++)) ; do
    ((res0 = res - 1))
    ((res1 = res + 1))
    if [[ $res -ne 1 ]] ; then
      echo "dihedral out dihedral.dat :$res0@C :$res@N :$res@CA :$res@C phi:$res" >> dih.in
    fi
    echo "dihedral out dihedral.dat :$res@N :$res@CA :$res@C  :$res1@N psi:$res" >> dih.in
  done
  # Test custom dihedral definition, ch1 2 4 5 6 8 9 11 12
  for res in 2 4 5 6 8 9 11 12 ; do
    echo "dihedral out dihedral2.dat :$res@N :$res@CA :$res@CB :$res@CG chi1:$res" >> dih.in
  done
  echo "multidihedral dihtype chi1:N:CA:CB:CG out custom.dat" >> dih.in
  echo "multidihedral out all.dat resrange 1-3" >> dih.in
  RunCpptraj "$UNITNAME"
  DoTest dihedral.dat multidih.dat
  DoTest dihedral2.dat custom.dat
  DoTest all.dat.save all.dat
fi

# Test nucleic acid chi
UNITNAME='Multidihedral Nucleotide CHI test'
CheckFor maxthreads 3
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ../adh026.3.pdb 
trajin ../adh026.3.pdb 
multidihedral chin out chin.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest chin.dat.save chin.dat
fi

# Test protein dihedrals 
UNITNAME='Multidihedral protein dihedrals test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm ARG.mol2
trajin ARG.mol2
multidihedral out arg.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest arg.dat.save arg.dat
fi

UNITNAME='Multidihedral cyclic molecule test'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > dih.in <<EOF
parm cyclic.mol2
trajin cyclic.mol2
multidihedral phi psi omega out cyclic.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest cyclic.dat.save cyclic.dat
fi

EndTest

exit 0
