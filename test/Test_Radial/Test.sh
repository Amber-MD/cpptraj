#!/bin/bash

. ../MasterTest.sh

CleanFiles radial.in Radial.agr cRadial.agr WatO-Trp4.agr WatO-Trp4.raw.agr \
           WatO-Trp4.byres.agr WatO-Trp.agr WatO-Trp.volume.agr \
           WatO-Glu5CD.agr noimage.WatO-Glu5CD.agr point.dat \
           point?.agr wat.origin.agr \
           watO-protein.agr watO-protein.raw.agr

TESTNAME='Radial tests'
Requires netcdf maxthreads 10

INPUT="-i radial.in"

UNITNAME='Radial test, non-orthogonal imaging'
cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2

radial WATO out watO-protein.agr rawrdf watO-protein.raw.agr :WAT@O ^1 0.5 10.0

EOF
RunCpptraj "$UNITNAME"
DoTest watO-protein.agr.save watO-protein.agr
DoTest watO-protein.raw.agr.save watO-protein.raw.agr

EndTest

exit 0
