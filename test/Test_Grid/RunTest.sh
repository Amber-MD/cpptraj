#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in out.dipole out.xplor out.dx out.dx.2 test.dx box.dx mask.dx \
           nonortho.dx triclinic.nc nonortho.pdb

CheckNetcdf
INPUT="ptraj.in"

# dipole
Dipole() {
  TOP="../tz2.ortho.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
rms first :1-13
center :1-13 mass origin 
image origin center familiar
dipole out.dipole 20 0.5 20 0.5 20 0.5 :WAT
EOF
  RunCpptraj "Dipole test"
  DoTest out.dipole.save out.dipole
  CheckTest
}

# grid
Grid() {
  TOP="../tz2.truncoct.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc
rms first :1-13
center :1-13 mass origin 
image origin center familiar
grid out.xplor 20 0.5 20 0.5 20 0.5 :WAT@O  
grid out.dx 20 0.5 20 0.5 20 0.5 :WAT@O
EOF
  RunCpptraj "Grid test"
  DoTest out.xplor.save out.xplor
  DoTest out.dx.save out.dx
  CheckTest
}

# grid dx read
GridDxRead() {
  TOP="../tz2.truncoct.parm7"
  cat > ptraj.in <<EOF
readdata out.dx.save
trajin ../tz2.truncoct.nc 1 1
grid out.dx.2 data out.dx.save :1588 opendx  
EOF
  RunCpptraj "OpenDX Grid read test"
  DoTest out.dx.save out.dx.2
  CheckTest
}

# Specified center
SpecifiedCenter() {
  TOP="../tz2.ortho.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc 1 10
autoimage origin
grid test.dx 34 .5 44 .5 36 .5 gridcenter 1.5 1.0 0.0 :WAT
EOF
  RunCpptraj "Grid with specified center test."
  DoTest test.dx.save test.dx
}

# Box center offset
BoxCenterOffset() {
  TOP="../tz2.ortho.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc 1 10
grid box.dx 30 .5 30 .5 30 .5 box :WAT
EOF
  RunCpptraj "Grid with box center offset test."
  DoTest box.dx.save box.dx
}

# Mask center offset
MaskCenterOffset() {
  TOP="../tz2.ortho.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc 1 10
grid mask.dx 30 .5 30 .5 30 .5 center :1-12@CA :WAT
EOF
  RunCpptraj "Grid with mask center offset test."
  DoTest mask.dx.save mask.dx
}

# Non-orthogonal grid
NonorthogonalGrid() {
  TOP="../tz2.truncoct.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc
reference ../tz2.truncoct.nc [REF]
autoimage triclinic
grid nonortho.dx boxref [REF] 50 50 50 :WAT@O pdb nonortho.pdb
#trajout triclinic.nc
EOF
  RunCpptraj "Non-orthogonal grid test."
  DoTest nonortho.dx.save nonortho.dx
}

Dipole
Grid
GridDxRead
SpecifiedCenter
BoxCenterOffset
MaskCenterOffset
NonorthogonalGrid

EndTest
exit 0
