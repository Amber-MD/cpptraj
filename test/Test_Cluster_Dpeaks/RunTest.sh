#!/bin/bash

. ../MasterTest.sh

CleanFiles dpeaks.in dpeaks.dat summary.dat info.dat

INPUT="-i dpeaks.in"
TESTNAME='Density peaks clustering test'
Requires netcdf
cat > dpeaks.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc
createcrd MyCrd
run

# Pass one - create dpeaks.dat
runanalysis cluster crdset MyCrd DPEAKS @CA dpeaks dvdfile dpeaks.dat epsilon 1.7 \
                    pairdist MyPairDist

# Manual
runanalysis cluster crdset MyCrd C0     @CA dpeaks dvdfile dpeaks.dat epsilon 1.7 \
                    pairdist MyPairDist \
                    choosepoints manual distancecut 2.4 densitycut 5.0 \
                    summary summary.dat info info.dat
#  dvdfile dpeaks.dat choosepoints manual distancecut 2.4 densitycut 0.0
EOF
RunCpptraj "$TESTNAME"
DoTest dpeaks.dat.save dpeaks.dat
DoTest summary.dat.save summary.dat
DoTest info.dat.save info.dat

EndTest
exit 0
