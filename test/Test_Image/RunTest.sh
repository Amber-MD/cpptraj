#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles image.in ortho.dat nonortho.dat image.crd image2.crd image3.crd image4.crd

INPUT="-i image.in"
TESTNAME='Imaging tests'
Requires netcdf maxthreads 10
# Test 1 - orthorhombic imaged distance
cat > image.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
distance image :507@O :499@O out ortho.dat
distance noimage :507@O :499@O out ortho.dat noimage
go
EOF
RunCpptraj "Orthorhombic imaged distance test."
DoTest ortho.dat.save ortho.dat

# Test 2 - nonorthorhombic imaged distance
cat > image.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 
distance image :1507@O :1676@O out nonortho.dat
distance noimage :1507@O :1676@O out nonortho.dat noimage
go
EOF
RunCpptraj "Non-orthorhombic imaged distance test."
DoTest nonortho.dat.save nonortho.dat

UNITNAME='Orthorhombic coordinate imaging test'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > image.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 2
center :2 origin
image origin center
trajout image.crd 
go
EOF
  RunCpptraj "$UNITNAME"
  DoTest image.crd.save image.crd
fi

UNITNAME='Nonorthorhombic coordinate imaging test'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > image.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
center :2-11
image center triclinic
trajout image2.crd
go
EOF
  RunCpptraj "$UNITNAME"
  DoTest image2.crd.save image2.crd
fi

UNITNAME='Nonorthorhombic coordinate imaging test with familiar'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > image.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
center :2-11 origin
image origin center familiar
trajout image3.crd
go
EOF
  RunCpptraj "$UNITNAME"
  DoTest image3.crd.save image3.crd
fi

UNITNAME='Nonorthorhombic coordinate imaging test with familiar and COM'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > image.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
center :2-11
image center familiar com :6
trajout image4.crd
go
EOF
  RunCpptraj "$UNITNAME"
  DoTest image4.crd.save image4.crd
fi
EndTest

exit 0
