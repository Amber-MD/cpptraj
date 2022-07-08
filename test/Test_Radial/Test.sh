#!/bin/bash

. ../MasterTest.sh

CleanFiles radial.in Radial.agr cRadial.agr WatO-Trp4.agr WatO-Trp4.raw.agr \
           WatO-Trp4.byres.agr WatO-Trp.agr WatO-Trp.volume.agr \
           WatO-Glu5CD.agr noimage.WatO-Glu5CD.agr point.dat \
           point?.agr wat.origin.agr \
           watO-protein.agr watO-protein.raw.agr \
           watO.agr watO.raw.agr

TESTNAME='Radial tests'
Requires netcdf maxthreads 10

INPUT="-i radial.in"

Radial_WatProt_Nonortho() {
  UNITNAME='Radial test (water-protein), non-orthogonal imaging'
  cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2

radial WATO out watO-protein.agr rawrdf watO-protein.raw.agr :WAT@O ^1 0.5 10.0

EOF
  RunCpptraj "$UNITNAME"
  DoTest watO-protein.agr.save watO-protein.agr
  DoTest watO-protein.raw.agr.save watO-protein.raw.agr
}

Radial_Wat_Nonortho() {
  UNITNAME='Radial test (water), non-orthogonal imaging'
  cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2

radial WATO out watO.agr rawrdf watO.raw.agr :WAT@O 0.5 10.0

EOF
  RunCpptraj "$UNITNAME"
  DoTest watO.agr.save watO.agr
  DoTest watO.raw.agr.save watO.raw.agr
}

Radial_WatProt_Nonortho
Radial_Wat_Nonortho

EndTest

exit 0
