#!/bin/bash

. ../MasterTest.sh

CleanFiles radial.in Radial.agr cRadial.agr WatO-Trp4.agr WatO-Trp4.raw.agr \
           WatO-Trp4.byres.agr WatO-Trp.agr WatO-Trp.volume.agr \
           WatO-Glu5CD.agr noimage.WatO-Glu5CD.agr point.dat \
           point?.agr wat.origin.agr tz2.WatO-Prot.agr tz2.WatO.agr

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

UNITNAME='Radial test, specified point'
CheckFor maxthreads 1
if [ $? -eq 0 ] ; then
  cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 1

#radial point0.agr    0.5 10.0 :WAT@O :5@CD
radial point1.agr    0.5 10.0 :WAT@O toxyz 11.5742,4.3807,-13.7675
#vector center :5@CD out point.dat
EOF
  RunCpptraj "$UNITNAME"
  DoTest point1.agr.save point1.agr
fi

UNITNAME='Radial test, specified point after centering'
cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc

autoimage origin
radial out wat.origin.agr 0.25 15.0 :WAT@O toxyz 0,0,0
EOF
RunCpptraj "$UNITNAME"
DoTest wat.origin.agr.save wat.origin.agr

UNITNAME='Radial test, water -> protein, non-orthogonal cell'
cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc
radial WatO-Prot :WAT@O ^1 out tz2.WatO-Prot.agr 0.25 15.0
EOF
RunCpptraj "$UNITNAME"
DoTest tz2.WatO-Prot.agr.save tz2.WatO-Prot.agr  

UNITNAME='Radial test, water -> water, non-orthogonal cell'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > radial.in <<EOF
parm ../tz2.truncoct.parm7
trajin ../tz2.truncoct.nc 1 2
radial WatO      :WAT@O    out tz2.WatO.agr      0.25 15.0
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2.WatO.agr.save tz2.WatO.agr
fi

UNITNAME='Radial test, water -> protein, orthogonal cell'
cat > radial.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
radial WatO-Prot :WAT@O ^1 out tz2ortho.WatO-Prot.agr 0.25 15.0
EOF
RunCpptraj "$UNITNAME"
DoTest tz2ortho.WatO-Prot.agr.save tz2ortho.WatO-Prot.agr  

UNITNAME='Radial test, water -> water, orthogonal cell'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > radial.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 2
radial WatO      :WAT@O    out tz2ortho.WatO.agr      0.25 15.0
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2ortho.WatO.agr.save tz2ortho.WatO.agr
fi

UNITNAME='Radial test, water -> protein, no imaging'
cat > radial.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc
radial WatO-Prot :WAT@O ^1 out tz2noimage.WatO-Prot.agr 0.25 15.0 noimage
EOF
RunCpptraj "$UNITNAME"
DoTest tz2noimage.WatO-Prot.agr.save tz2noimage.WatO-Prot.agr  

UNITNAME='Radial test, water -> water, no imaging'
CheckFor maxthreads 2
if [ $? -eq 0 ] ; then
  cat > radial.in <<EOF
parm ../tz2.ortho.parm7
trajin ../tz2.ortho.nc 1 2
radial WatO      :WAT@O    out tz2noimage.WatO.agr      0.25 15.0 noimage
EOF
  RunCpptraj "$UNITNAME"
  DoTest tz2noimage.WatO.agr.save tz2noimage.WatO.agr
fi

EndTest

exit 0
