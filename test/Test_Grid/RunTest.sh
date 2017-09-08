#!/bin/bash

. ../MasterTest.sh

CleanFiles ptraj.in out.dipole out.xplor out.dx out.dx.2 test.dx box.dx mask.dx \
           nonortho.dx triclinic.nc nonortho.pdb bounds.dat bounds.mol2 \
           bounds.xplor avg.mol2 nonortho_wrap.dx

TESTNAME='Grid tests'
Requires netcdf maxthreads 10
INPUT="ptraj.in"

# dipole
Dipole() {
  TOP="../tz2.ortho.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
autoimage origin
rms first :1-13
dipole out.dipole 20 0.5 20 0.5 20 0.5 :WAT
EOF
  RunCpptraj "Dipole test"
  DoTest out.dipole.save out.dipole
}

# grid
Grid() {
  TOP="../tz2.truncoct.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc
autoimage origin
rms first :1-13
average avg.mol2 :1-13 
grid out.xplor 20 0.5 20 0.5 20 0.5 :WAT@O name XPLOR
grid out.dx 20 0.5 20 0.5 20 0.5 :WAT@O
EOF
  RunCpptraj "Grid test"
  DoTest out.xplor.save out.xplor
  DoTest out.dx.save out.dx
}

# grid dx read
GridDxRead() {
  TOP="../tz2.truncoct.parm7"
  cat > ptraj.in <<EOF
readdata out.dx.save
trajin ../tz2.truncoct.nc
autoimage origin
rms first :1-13
grid out.dx.2 data out.dx.save @CA opendx  
EOF
  RunCpptraj "OpenDX Grid read test"
  DoTest out.dx.2.save out.dx.2
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

# Non-orthogonal grid, points centered on bins
NonorthoGridBinCenter() {
  TOP="../tz2.truncoct.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.truncoct.nc
reference ../tz2.truncoct.nc [REF]
autoimage triclinic
grid nonortho_wrap.dx boxref [REF] 20 20 20 :WAT@O gridwrap name BC
EOF
  RunCpptraj "Non-orthogonal grid centered on bin center and wrapped test"
  DoTest nonortho_wrap.dx.save nonortho_wrap.dx
}

# Generate grid from bounds
Bounds() {
  TOP="../tz2.ortho.parm7"
  cat > ptraj.in <<EOF
trajin ../tz2.ortho.nc
autoimage
rms first :1-13&!@H= mass
bounds :1-13 dx .5 name MyGrid out bounds.dat
#average bounds.mol2 :1-13
createcrd MyCoords
run
crdaction MyCoords grid bounds.xplor data MyGrid :WAT@O
EOF
  RunCpptraj "Grid generation from 'bounds' test."
  DoTest bounds.dat.save bounds.dat
  DoTest bounds.xplor.save bounds.xplor 
}

Dipole
Grid
GridDxRead
SpecifiedCenter
BoxCenterOffset
MaskCenterOffset
NonorthogonalGrid
NonorthoGridBinCenter
Bounds

EndTest
exit 0
