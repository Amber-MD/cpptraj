#!/bin/bash

# Test of IRED 
. ../MasterTest.sh

CleanFiles orderparam ired.vec noe v0.cjt v0.cmt ired.in v0 plateau.norm.dat cjt?.dat \
           matrix_ds2.dat

TESTNAME='IRED vector/matrix test'
Requires mathlib

TOP="1IEE_A_prot.prmtop"
INPUT="ired.in"
cat > $INPUT <<EOF
trajin 1IEE_A_test.mdcrd
for atoms Natom inmask :2-129@N&!:PRO atoms Hatom inmask :2-129@H v=1;v++
  vector v\$v \$Natom ired \$Hatom
done

matrix ired name matired order 2
diagmatrix matired out ired.vec name ired.vec vecs \$v
ired order 2 modes ired.vec relax freq 500.0 NHdist 1.02 \
     tstep 1.0 tcorr 10000.0 norm name MyIred \
     out v0 noefile noe ds2matrix matrix_ds2.dat 
datafile v0.cmt noheader
datafile v0.cjt noheader
create plateau.norm.dat MyIred[Plateau] noxcol prec 12.8 noheader
create orderparam MyIred[S2] prec 10.5
run
writedata cjt1.dat MyIred[Cj(t)]:0 
EOF

RunCpptraj "$TESTNAME"
DoTest ired.vec.save ired.vec
DoTest v0.cmt.new.norm.save v0.cmt
DoTest plateau.norm.dat.save plateau.norm.dat
DoTest v0.cjt.new.save v0.cjt
DoTest orderparam.save orderparam
DoTest matrix_ds2.dat.save matrix_ds2.dat
DoTest noe.save noe

#DoTest v0.cjt.save v0.cjt
#DoTest v0.cmt.save v0.cmt
EndTest

exit 0
