#!/bin/bash

. ../MasterTest.sh

CleanFiles hist.in hist.gnu hist.agr freeE.gnu norm.gnu hist.dx 3D.dat \
           phipsihist.agr phipsikde.agr
INPUT="-i hist.in"

# 1D and 2D histogram tests
UNITNAME='1D/2D Histogram Analysis test'
CheckFor netcdf
if [ $? -eq 0 ] ; then
  cat > hist.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

distance d1 :1 :10 
distance d2 :10 :13  

dihedral phi6 :5@C :6@N :6@CA :6@C 
dihedral psi6 :6@N :6@CA :6@C :7@N 

hist d1 d2,8,10,0.1,* min 9.0 max 26.0 step 0.5 out hist.gnu
hist name D1hist d1,9,26,0.5 out hist.agr
hist phi6 psi6 min -180 max 180 bins 72 out freeE.gnu free 300.0
hist phi6 psi6 min -180 max 180 bins 72 out norm.gnu norm
multihist phi6 psi6 min -180 max 180 bins 72 out phipsihist.agr 
multihist phi6 psi6 min -180 max 180 bins 72 kde out phipsikde.agr 
EOF
  RunCpptraj "$UNITNAME"
  DoTest hist.gnu.save hist.gnu
  DoTest hist.agr.save hist.agr
  DoTest freeE.gnu.save freeE.gnu
  DoTest norm.gnu.save norm.gnu
  DoTest phipsihist.agr.save phipsihist.agr
  DoTest phipsikde.agr.save phipsikde.agr
fi

# 3D histogram test
cat > 3D.dat <<EOF
1.0 1.0 1.0
2.0 1.0 1.0
2.0 2.0 2.0
3.0 3.0 3.0
4.0 4.0 4.0
EOF
cat > hist.in <<EOF
readdata 3D.dat 
runanalysis hist 3D.dat:1 3D.dat:2 3D.dat:3 min 0.5 max 4.5 step 1.0 out hist.dx
EOF
RunCpptraj "3D histogram test."
DoTest hist.dx.save hist.dx

EndTest

exit 0
