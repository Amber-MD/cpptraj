#!/bin/bash

. ../MasterTest.sh

CleanFiles vector.in v6and7.dat rex-d.dat MD.ene.dat out.dx out.dat \
           truncoct.dat out?.dx append.dx append.dat temp.dat \
           truncsparse.dat temp2.dat truncoct.dat.save sparse.dat \
           xyz.dat

TESTNAME='Read data tests'

INPUT="-i vector.in"
# Test read/append of vector dataset
UNITNAME='Read vector data test'
cat > vector.in <<EOF
readdata ../Test_Vector/vtest.dat.6.save vector name v6and7
readdata ../Test_Vector/vtest.dat.7.save vector name v6and7
writedata v6and7.dat v6and7
EOF
RunCpptraj "$UNITNAME"
DoTest v6and7.dat.save v6and7.dat

UNITNAME='Read Amber output test'
cat > vector.in <<EOF
readdata md.initial.out md.restart.out name MD
writedata MD.ene.dat MD[*] prec 14.4
EOF
RunCpptraj "$UNITNAME"
DoTest MD.ene.dat.save MD.ene.dat

UNITNAME='Read CHARMM output test'
cat > vector.in <<EOF
readdata rex-d.out_0 name ENE
writedata rex-d.dat ENE[*]
EOF
RunCpptraj "$UNITNAME"
DoTest rex-d.dat.save rex-d.dat

# Make sure we can read data in and write it out properly
UNITNAME='3D data read/write test'
cat > vector.in <<EOF
readdata ../Test_Grid/out.dx.save name grid
writedata out.dx grid
EOF
RunCpptraj "$UNITNAME"
DoTest ../Test_Grid/out.dx.save out.dx

# Read data as DX, write standard
UNITNAME='Standard 3D data write test'
cat > vector.in <<EOF
readdata ../Test_Grid/out.dx.save name grid
list datasets
writedata out.dat grid
EOF
RunCpptraj "$UNITNAME"
# Read data standard, write DX 
UNITNAME='Standard 3D data read test'
cat > vector.in <<EOF
readdata out.dat read3d dims 20,20,20 origin -5,-5,-5 delta .5,.5,.5 name MyGrid prec dbl
list datasets
writedata out2.dx MyGrid
EOF
RunCpptraj "$UNITNAME"
DoTest ../Test_Grid/out.dx.save out2.dx

# Test appending grid data
UNITNAME='Standard 3D data read/append test'
cat > vector.in <<EOF
readdata ../Test_Grid/out.dx.save name grid
writedata out.dat grid
readdata out.dat read3d dims 20,20,20 origin -5,-5,-5 delta .5,.5,.5 name grid
writedata append.dat grid sparse
EOF
RunCpptraj "$UNITNAME"
DoTest append.dat.save append.dat

# Test write/read sparse grid
UNITNAME='Standard 3D data, sparse'
cat > vector.in <<EOF
readdata ../Test_Grid/out.dx.save name grid
writedata sparse.dat grid sparse
# Read the sparse set back in
readdata sparse.dat read3d dims 20,20,20 origin -5,-5,-5 delta .5,.5,.5 name grid2
writedata out3.dx grid2
EOF
RunCpptraj "$UNITNAME"
DoTest ../Test_Grid/out.dx.save out3.dx

# Nonorthogonal grid write
UNITNAME='Nonorthogonal grid standard write/read tests.'
CheckFor netcdf maxthreads 1
if [ $? -eq 0 ] ; then
  cat > vector.in <<EOF
parm ../tz2.truncoct.parm7
reference ../tz2.truncoct.nc [first]
trajin ../tz2.truncoct.nc 1 1
autoimage triclinic
grid boxref [first] 42 42 42 :WAT@O normdensity name MyGrid
run
# Write data as sparse
writedata truncsparse.dat MyGrid sparse
# Write data as complete
writedata truncoct.dat.save MyGrid
clear datasets
# Test that we can read complete nonortho data back in
readdata truncoct.dat.save read3d name MyGrid
writedata temp.dat MyGrid
clear datasets
# Test that we can read sparse nonortho data
readdata truncsparse.dat read3d name MyGrid
writedata temp2.dat MyGrid
EOF
  RunCpptraj "$UNITNAME"
  DoTest truncsparse.dat.save truncsparse.dat
  DoTest truncoct.dat.save temp.dat
  DoTest truncoct.dat.save temp2.dat
fi

# Read specific columns
UNITNAME='Read specific columns test.'
cat > vector.in <<EOF
readdata ../Test_Diffusion/diff.1.dat.save onlycols 1,6-9 index 1 name XYZ
writedata xyz.dat XYZ
EOF
RunCpptraj "$UNITNAME"
DoTest xyz.dat.save xyz.dat

EndTest
  
exit 0
