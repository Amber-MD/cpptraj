#!/bin/bash

. ../MasterTest.sh

CleanFiles radial.in Radial.agr cRadial.agr WatO-Trp4.agr WatO-Trp4.raw.agr
CheckNetcdf
cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10

radial Radial.agr    0.5 10.0 :5@CD  :WAT@O
radial cRadial.agr   0.5 10.0 :5     :WAT@O center1
radial WatO-Trp4.agr 0.5 10.0 :WAT@O :4&!@C,O,CA,HA,N,H center2 \
       intrdf WatO-Trp4.raw.agr rawrdf WatO-Trp4.raw.agr
EOF

INPUT="-i radial.in"
RunCpptraj "Radial Test"
DoTest Radial.agr.save Radial.agr
DoTest cRadial.agr.save cRadial.agr
DoTest WatO-Trp4.agr.save WatO-Trp4.agr
DoTest WatO-Trp4.raw.agr.save WatO-Trp4.raw.agr

EndTest

exit 0
