#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles image.in ortho.dat nonortho.dat image.crd image2.crd image3.crd image4.crd

# Test - orthorhombic imaged distance
cat > image.in <<EOF
noprogress
parm Ala10.parm7
trajin Ala10.nc
distance image @2 @107 out ortho.dat
distance noimage @2 @107 out ortho.dat noimage
go
EOF
INPUT="-i image.in"
RunCpptraj "Orthorhombic imaged distance test."

# Test - nonorthorhombic imaged distance
cat > image.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
distance image :49@HZ3 :258@HD22 out nonortho.dat
distance noimage :49@HZ3 :258@HD22 out nonortho.dat noimage
go
EOF
INPUT="-i image.in"
RunCpptraj "Non-orthorhombic imaged distance test."

DoTest ortho.dat.save ortho.dat
DoTest nonortho.dat.save nonortho.dat
CheckTest

# Test - Orthorhombic coordinate imaging 
cat > image.in <<EOF
noprogress
parm Ala10.parm7
trajin Ala10.nc
center :2
image origin center
trajout image.crd 
go
EOF
INPUT="-i image.in"
RunCpptraj "Orthorhombic coordinate imaging test."

# Test - Nonorthorhombic coordinate imaging
cat > image.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
center :190-211
image origin center
trajout image2.crd
go
EOF
INPUT="-i image.in"
RunCpptraj "Nonorthorhombic coordinate imaging test."

# Test - Nonorthorhombic coordinate imaging with familiar
cat > image.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
center :190-211
image origin center familiar
trajout image3.crd
go
EOF
INPUT="-i image.in"
RunCpptraj "Nonorthorhombic coordinate imaging test with familiar."

# Test - Nonorthorhombic coordinate imaging test with familiar and COM
cat > image.in <<EOF
noprogress
parm ../ChainA-tip3p.parm7
trajin ../run0.nc 1 10
center :190-211
image origin center familiar com :100
trajout image4.crd
go
EOF
INPUT="-i image.in"
RunCpptraj "Nonorthorhombic coordinate imaging test with familiar and COM."

DoTest image.crd.save image.crd
DoTest image2.crd.save image2.crd
DoTest image3.crd.save image3.crd
DoTest image4.crd.save image4.crd
CheckTest

EndTest

exit 0
