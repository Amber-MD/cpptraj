#!/bin/bash

. ../MasterTest.sh

CleanFiles tica.in ticadebug.dat M?.dat ticadebug.m?.dat \
           tz2.*.dat crd.*.dat dih.*.dat

INPUT='-i tica.in'

#MASS=mass
#cat > tica.in <<EOF
#parm ../tz2.parm7
#loadcrd ../tz2.crd name TZ2
#
#set MASK = :1-3@N,CA,C
#
#crdaction TZ2 matrix name M1 mwcovar out M1.dat \$MASK $MASS
#
#runanalysis tica crdset TZ2 mask \$MASK lag 1 debugc0 ticadebug.m1.dat $MASS
#
#crdaction TZ2 matrix name M2 mwcovar out M2.dat :1-3@N,CA :1-3@C,O $MASS
#
#runanalysis tica crdset TZ2 mask :1-3@N,CA mask2 :1-3@C,O lag 1 debugct ticadebug.m2.dat $MASS
#
#list dataset
#EOF
cat > tica.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
#createcrd TZ2
distance d1.12 :1@CA :12@CA
distance d2.11 :2@CA :11@CA
distance d3.10 :3@CA :10@CA
distance d4.9  :4@CA :9@CA
distance d5.8  :5@CA :8@CA
rms R0 first @CA
tica data R0 data d* name TICA lag 10 out tz2.ticamodes.dat cumvarout tz2.cumvar.dat
tica data R0 data d* name TICA2 lag 10 map commute out tz2.tica2modes.dat
run
writedata tz2.raw.dat R0 d*

#projection Evec out tz2.project.dat evecs TICA data R0 data d*
runanalysis projectdata name Evec out tz2.project.dat evecs TICA data R0 data d*
run
EOF
RunCpptraj "TICA test, 1D data sets."
DoTest tz2.ticamodes.dat.save tz2.ticamodes.dat
DoTest tz2.cumvar.dat.save tz2.cumvar.dat
DoTest tz2.project.dat.save tz2.project.dat
DoTest tz2.tica2modes.dat.save tz2.tica2modes.dat

cat > tica.in <<EOF
parm ../tz2.parm7
loadcrd ../tz2.crd name MyCrd
tica crdset MyCrd mask @CA name TICA lag 10 out crd.ticamodes.dat cumvarout crd.cumvar.dat
run
# DEBUG
#crdaction MyCrd strip !@CA
#runanalysis tica crdset MyCrd name TICA2 lag 10 out crd.tica2modes.dat cumvarout crd.cumvar2.dat
EOF
RunCpptraj "TICA test, coordinates."
DoTest crd.ticamodes.dat.save crd.ticamodes.dat
DoTest crd.cumvar.dat.save crd.cumvar.dat

cat > tica.in <<EOF
parm ../tz2.parm7
trajin ../tz2.crd
dihedral d1 :1@C :2@N :2@CA :2@C
dihedral d2 :2@N :2@CA :2@C :3@N
createcrd MyCrd
tica data d1 data d2 name TICA lag 10 out dih.ticamodes.dat cumvarout dih.cumvar.dat
run
EOF
RunCpptraj "TICA test, dihedrals."
DoTest dih.ticamodes.dat.save dih.ticamodes.dat
DoTest dih.cumvar.dat.save dih.cumvar.dat

#diff M1.dat ticadebug.m1.dat
#diff M2.dat ticadebug.m2.dat

EndTest
exit 0
