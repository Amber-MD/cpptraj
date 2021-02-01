#!/bin/bash

. ../MasterTest.sh

CleanFiles radial.in Radial.agr cRadial.agr WatO-Trp4.agr WatO-Trp4.raw.agr \
           WatO-Trp4.byres.agr WatO-Trp.agr WatO-Trp.volume.agr \
           WatO-Glu5CD.agr noimage.WatO-Glu5CD.agr

TESTNAME='Radial tests'
Requires netcdf maxthreads 10

INPUT="-i radial.in"

UNITNAME='Radial test, non-orthogonal imaging'
cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 10

radial Radial.agr    0.5 10.0 :5@CD  :WAT@O
radial cRadial.agr   0.5 10.0 :5     :WAT@O center1
radial WatO-Trp4.agr 0.5 10.0 :WAT@O :4&!@C,O,CA,HA,N,H center2 \
       intrdf WatO-Trp4.raw.agr rawrdf WatO-Trp4.raw.agr
radial WatO-Trp4.byres.agr 0.5 10.0 :WAT@O :4&!@C,O,CA,HA,N,H byres2
radial out WatO-Trp.agr 0.5 10.0 :WAT@O :TRP byres2
radial out WatO-Trp.volume.agr 0.5 10.0 :WAT@O :TRP volume  
EOF
RunCpptraj "$UNITNAME"
DoTest Radial.agr.save Radial.agr
DoTest cRadial.agr.save cRadial.agr
DoTest WatO-Trp4.agr.save WatO-Trp4.agr
DoTest WatO-Trp4.raw.agr.save WatO-Trp4.raw.agr
DoTest WatO-Trp4.agr.save WatO-Trp4.byres.agr
DoTest WatO-Trp.agr.save WatO-Trp.agr
DoTest WatO-Trp.volume.agr.save WatO-Trp.volume.agr

UNITNAME='Radial test, orthogonal and no imaging'
cat > radial.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
radial O_CD :WAT@O :5@CD out WatO-Glu5CD.agr 0.5 20.0
radial O_CD.noimage :WAT@O :5@CD out noimage.WatO-Glu5CD.agr 0.5 20.0 noimage
EOF
RunCpptraj "$UNITNAME"
DoTest WatO-Glu5CD.agr.save WatO-Glu5CD.agr
DoTest noimage.WatO-Glu5CD.agr.save noimage.WatO-Glu5CD.agr

EndTest

exit 0
