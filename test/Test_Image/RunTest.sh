#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles image.in ortho.dat nonortho.dat image.crd image2.crd image3.crd image4.crd

INPUT="-i image.in"
RequiresNetcdf "Imaging tests"
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

MaxThreads 2 "Coordinate imaging tests"
if [ $? -ne 0 ] ; then
  SkipCheck "Coordinate imaging tests"
else
  # Test 3 - Orthorhombic coordinate imaging 
  cat > image.in <<EOF
noprogress
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 2
center :2 origin
image origin center
trajout image.crd 
go
EOF
  RunCpptraj "Orthorhombic coordinate imaging test."
  DoTest image.crd.save image.crd

  # Test 4 - Nonorthorhombic coordinate imaging
  cat > image.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
center :2-11
image center triclinic
trajout image2.crd
go
EOF
  RunCpptraj "Nonorthorhombic coordinate imaging test."
  DoTest image2.crd.save image2.crd

  # Test - Nonorthorhombic coordinate imaging with familiar
  cat > image.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
center :2-11 origin
image origin center familiar
trajout image3.crd
go
EOF
  RunCpptraj "Nonorthorhombic coordinate imaging test with familiar."
  DoTest image3.crd.save image3.crd

  # Test - Nonorthorhombic coordinate imaging test with familiar and COM
cat > image.in <<EOF
noprogress
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
center :2-11
image center familiar com :6
trajout image4.crd
go
EOF
  RunCpptraj "Nonorthorhombic coordinate imaging test with familiar and COM."
  DoTest image4.crd.save image4.crd
fi
EndTest

exit 0
