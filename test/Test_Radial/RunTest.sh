#!/bin/bash

. ../MasterTest.sh

CleanFiles radial.in Radial.agr cRadial.agr
CheckNetcdf
cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10

radial Radial.agr 0.5 10.0 :5@CD :WAT@O 
radial cRadial.agr 0.5 10.0 :5 :WAT@O center1
EOF

INPUT="-i radial.in"
RunCpptraj "Radial Test"
DoTest Radial.agr.save Radial.agr
DoTest cRadial.agr.save cRadial.agr
CheckTest

EndTest

exit 0
